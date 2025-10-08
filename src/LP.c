// GLPK-free dense primal simplex (two-phase), with per-ballot rescale of W if totals mismatch.

#include "globals.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ---- BLAS helpers ----
static inline void row_scale(int n, double *row, double alpha)
{
    int inc = 1;
    F77_CALL(dscal)(&n, &alpha, row, &inc);
}
static inline void row_axpy(int n, double a, const double *x, double *y)
{
    int inc = 1;
    F77_CALL(daxpy)(&n, &a, x, &inc, y, &inc);
}
static inline void row_copy(int n, const double *src, double *dst)
{
    int inc = 1;
    F77_CALL(dcopy)(&n, src, &inc, dst, &inc);
}
static void *xcalloc(size_t n, size_t sz)
{
    void *p = calloc(n, sz);
    if (!p)
        error("Allocation error");
    return p;
}
static inline int idxRC(int r, int ld, int c)
{
    return r * ld + c;
}

// ---- Simplex tab ----
typedef struct
{
    int m, n, ncols;
    double *T;
    int *basis;
    double tol;
} SimplexTab;

static void pivot_at(SimplexTab *tab, int leave, int enter)
{
    int m = tab->m, ncols = tab->ncols;
    double *prow = &tab->T[idxRC(leave, ncols, 0)];
    double piv = tab->T[idxRC(leave, ncols, enter)];
    if (fabs(piv) < 1e-14)
        return;
    row_scale(ncols, prow, 1.0 / piv);
    for (int i = 0; i <= m; ++i)
    {
        if (i == leave)
            continue;
        double *row = &tab->T[idxRC(i, ncols, 0)];
        double a = row[enter];
        if (fabs(a) > 0)
            row_axpy(ncols, -a, prow, row);
    }
    tab->basis[leave] = enter;
}
static int choose_entering_bland(const SimplexTab *t)
{
    const double *obj = &t->T[idxRC(t->m, t->ncols, 0)];
    for (int j = 0; j < t->n; ++j)
        if (obj[j] < -t->tol)
            return j;
    return -1;
}
static int choose_leaving(const SimplexTab *t, int enter)
{
    int m = t->m, ncols = t->ncols;
    const double *T = t->T;
    int leave = -1;
    double best = HUGE_VAL;
    for (int i = 0; i < m; ++i)
    {
        double a = T[idxRC(i, ncols, enter)];
        if (a > t->tol)
        {
            double rhs = T[idxRC(i, ncols, t->n)];
            double ratio = rhs / a;
            if (ratio < best - 1e-12 ||
                (fabs(ratio - best) <= 1e-12 && t->basis[i] > (leave >= 0 ? t->basis[leave] : -1)))
            {
                best = ratio;
                leave = i;
            }
        }
    }
    return leave;
}
static int simplex_solve_tab(SimplexTab *t)
{
    while (1)
    {
        int e = choose_entering_bland(t);
        if (e < 0)
            return 0;
        int l = choose_leaving(t, e);
        if (l < 0)
            return 1;
        pivot_at(t, l, e);
    }
}

// ---- Phase I/II ----
static void build_phase1_full(const double *A, const double *b, int m, int n, SimplexTab *t)
{
    int n_art = m;
    t->m = m;
    t->n = n + n_art;
    t->ncols = t->n + 1;
    t->T = (double *)xcalloc((size_t)(m + 1) * (t->ncols), sizeof(double));
    t->basis = (int *)xcalloc(m, sizeof(int));
    t->tol = 1e-10;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
            t->T[idxRC(i, t->ncols, j)] = A[idxRC(i, n, j)];
        int ja = n + i;
        t->T[idxRC(i, t->ncols, ja)] = 1.0;
        t->T[idxRC(i, t->ncols, t->n)] = b[i];
        t->basis[i] = ja;
    }
    double *obj = &t->T[idxRC(m, t->ncols, 0)];
    for (int j = n; j < n + n_art; ++j)
        obj[j] = 1.0;
    for (int i = 0; i < m; ++i)
    {
        double *row = &t->T[idxRC(i, t->ncols, 0)];
        row_axpy(t->ncols, -1.0, row, obj);
    }
}
static double phase1_artificial_value(const SimplexTab *p)
{
    int n_orig = p->n - p->m, ncols = p->ncols;
    double s = 0;
    for (int i = 0; i < p->m; ++i)
    {
        int b = p->basis[i];
        if (b >= n_orig)
            s += p->T[idxRC(i, ncols, p->n)];
    }
    return s;
}
static int phase1_eliminate_artificials(SimplexTab *p)
{
    int m = p->m, n_orig = p->n - m, ncols = p->ncols;
    for (int i = 0; i < m; ++i)
    {
        int b = p->basis[i];
        if (b >= n_orig)
        {
            int enter = -1;
            double best = 0;
            for (int j = 0; j < n_orig; ++j)
            {
                double a = p->T[idxRC(i, ncols, j)];
                if (fabs(a) > best + 1e-12)
                {
                    best = fabs(a);
                    enter = j;
                }
            }
            if (enter >= 0)
                pivot_at(p, i, enter);
            else
            {
                double rhs = p->T[idxRC(i, ncols, p->n)];
                if (fabs(rhs) > 1e-9)
                    return 1;
            }
        }
    }
    return 0;
}
static void set_phase2_objective(SimplexTab *t, const double *c, int n_orig)
{
    double *obj = &t->T[idxRC(t->m, t->ncols, 0)];
    for (int j = 0; j < t->n; ++j)
        obj[j] = (j < n_orig ? c[j] : 0);
    for (int i = 0; i < t->m; ++i)
    {
        int b = t->basis[i];
        if (b >= 0)
        {
            double coef = obj[b];
            if (fabs(coef) > 0)
            {
                double *row = &t->T[idxRC(i, t->ncols, 0)];
                row_axpy(t->ncols, -coef, row, obj);
            }
        }
    }
}
static void build_phase2_from_phase1(SimplexTab *p1, const double *c, SimplexTab *p2)
{
    int m = p1->m, n_orig = p1->n - m;
    if (phase1_eliminate_artificials(p1))
        p1->T[idxRC(m, p1->ncols, p1->n)] = HUGE_VAL;
    set_phase2_objective(p1, c, n_orig);
    p2->m = m;
    p2->n = n_orig;
    p2->ncols = n_orig + 1;
    p2->T = (double *)xcalloc((size_t)(m + 1) * (p2->ncols), sizeof(double));
    p2->basis = (int *)xcalloc(m, sizeof(int));
    p2->tol = 1e-10;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n_orig; ++j)
            p2->T[idxRC(i, p2->ncols, j)] = p1->T[idxRC(i, p1->ncols, j)];
        p2->T[idxRC(i, p2->ncols, p2->n)] = p1->T[idxRC(i, p1->ncols, p1->n)];
    }
    for (int j = 0; j < n_orig; ++j)
        p2->T[idxRC(m, p2->ncols, j)] = p1->T[idxRC(m, p1->ncols, j)];
    for (int i = 0; i < m; ++i)
    {
        int b = p1->basis[i];
        p2->basis[i] = (b >= 0 && b < n_orig) ? b : -1;
        if (p2->basis[i] >= 0)
        {
            double piv = p2->T[idxRC(i, p2->ncols, p2->basis[i])];
            if (fabs(piv) > 1e-12)
                pivot_at(p2, i, p2->basis[i]);
        }
    }
}
static void extract_solution(const SimplexTab *p2, double *x)
{
    int m = p2->m, n = p2->n, ncols = p2->ncols;
    memset(x, 0, sizeof(double) * n);
    for (int i = 0; i < m; ++i)
    {
        int b = p2->basis[i];
        if (0 <= b && b < n)
            x[b] = p2->T[idxRC(i, ncols, n)];
    }
}
static int simplex_solve_dense_equalities(double *A, double *b, const double *c, int m, int n, double *x)
{
    for (int i = 0; i < m; ++i)
    {
        if (b[i] < 0)
        {
            row_scale(n, &A[idxRC(i, n, 0)], -1.0);
            b[i] = -b[i];
        }
    }
    SimplexTab P1 = {0};
    build_phase1_full(A, b, m, n, &P1);
    int rc = simplex_solve_tab(&P1);
    double art = phase1_artificial_value(&P1);
    if (rc == 1 || art > 1e-7)
    {
        free(P1.T);
        free(P1.basis);
        return 1;
    }
    SimplexTab P2 = {0};
    build_phase2_from_phase1(&P1, c, &P2);
    free(P1.T);
    free(P1.basis);
    rc = simplex_solve_tab(&P2);
    if (rc == 1)
    {
        free(P2.T);
        free(P2.basis);
        return 2;
    }
    extract_solution(&P2, x);
    free(P2.T);
    free(P2.basis);
    return 0;
}

// ---- Layout ----
typedef struct
{
    int G, C, nQ, nY, nS, nVars, nAbs, nGrp, nCat, nRows, offQ, offY, offS;
} Layout;
static Layout mk_layout(int G, int C)
{
    Layout L = {.G = G, .C = C};
    L.nQ = G * C;
    L.nY = G * C;
    L.nAbs = 2 * G * C;
    L.nS = L.nAbs;
    L.nGrp = G;
    L.nCat = C;
    L.nVars = L.nQ + L.nY + L.nS;
    L.nRows = L.nAbs + L.nGrp + L.nCat;
    L.offQ = 0;
    L.offY = L.offQ + L.nQ;
    L.offS = L.offY + L.nY;
    return L;
}
static inline int col_q(const Layout *L, int g, int c)
{
    return L->offQ + g * L->C + c;
}
static inline int col_y(const Layout *L, int g, int c)
{
    return L->offY + g * L->C + c;
}
static inline int col_s1(const Layout *L, int g, int c)
{
    return L->offS + (g * L->C + c) * 2 + 0;
}
static inline int col_s2(const Layout *L, int g, int c)
{
    return L->offS + (g * L->C + c) * 2 + 1;
}
static inline int row_abs1(const Layout *L, int g, int c)
{
    return (g * L->C + c) * 2 + 0;
}
static inline int row_abs2(const Layout *L, int g, int c)
{
    return (g * L->C + c) * 2 + 1;
}
static inline int row_grp(const Layout *L, int g)
{
    return L->nAbs + g;
}
static inline int row_cand(const Layout *L, int c)
{
    return L->nAbs + L->nGrp + c;
}

// ---- Solve per ballot ----
static int solve_ballot_simplex(EMContext *ctx, int b, int mode)
{
    int G = ctx->G, C = ctx->C;
    Layout L = mk_layout(G, C);
    double sum_w = 0, sum_x = 0;
    for (int g = 0; g < G; ++g)
        sum_w += MATRIX_AT(ctx->W, b, g);
    for (int c = 0; c < C; ++c)
        sum_x += MATRIX_AT(ctx->X, c, b);
    double rel = (fmax(fabs(sum_w), fabs(sum_x)) > 0) ? fabs(sum_w - sum_x) / fmax(fabs(sum_w), fabs(sum_x)) : 0;
    double alpha = (rel > 1e-12 && sum_w > 0) ? (sum_x / sum_w) : 1.0;
    double *Weff = (double *)xcalloc(G, sizeof(double));
    for (int g = 0; g < G; ++g)
        Weff[g] = MATRIX_AT(ctx->W, b, g) * alpha;

    double *A = (double *)xcalloc((size_t)L.nRows * (size_t)L.nVars, sizeof(double));
    double *rhs = (double *)xcalloc(L.nRows, sizeof(double));
    double *obj = (double *)xcalloc(L.nVars, sizeof(double));
    for (int g = 0; g < G; ++g)
        for (int c = 0; c < C; ++c)
            obj[col_y(&L, g, c)] = 1.0;

    for (int g = 0; g < G; ++g)
    {
        double w = Weff[g];
        for (int c = 0; c < C; ++c)
        {
            double t = Q_3D(ctx->q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES);
            double coef = (mode ? w : 1.0);
            int r1 = row_abs1(&L, g, c), r2 = row_abs2(&L, g, c);
            A[idxRC(r1, L.nVars, col_y(&L, g, c))] = 1.0;
            A[idxRC(r1, L.nVars, col_q(&L, g, c))] = -coef;
            A[idxRC(r1, L.nVars, col_s1(&L, g, c))] = -1.0;
            rhs[r1] = -coef * t;
            A[idxRC(r2, L.nVars, col_y(&L, g, c))] = 1.0;
            A[idxRC(r2, L.nVars, col_q(&L, g, c))] = coef;
            A[idxRC(r2, L.nVars, col_s2(&L, g, c))] = -1.0;
            rhs[r2] = coef * t;
        }
    }
    for (int g = 0; g < G; ++g)
    {
        int r = row_grp(&L, g);
        for (int c = 0; c < C; ++c)
            A[idxRC(r, L.nVars, col_q(&L, g, c))] = 1.0;
        rhs[r] = 1.0;
    }
    for (int c = 0; c < C; ++c)
    {
        int r = row_cand(&L, c);
        double xcb = MATRIX_AT(ctx->X, c, b);
        for (int g = 0; g < G; ++g)
            A[idxRC(r, L.nVars, col_q(&L, g, c))] = Weff[g];
        rhs[r] = xcb;
    }

    double *x = (double *)xcalloc(L.nVars, sizeof(double));
    int code = simplex_solve_dense_equalities(A, rhs, obj, L.nRows, L.nVars, x);
    if (code != 0)
    {
        for (int g = 0; g < G; ++g)
            for (int c = 0; c < C; ++c)
                Q_3D(ctx->q, b, g, c, G, C) = 0.0;
    }
    else
    {
        for (int g = 0; g < G; ++g)
            for (int c = 0; c < C; ++c)
            {
                double q = x[col_q(&L, g, c)];
                if (q < 0 && q > -1e-12)
                    q = 0;
                Q_3D(ctx->q, b, g, c, G, C) = q;
            }
    }
    free(Weff);
    free(A);
    free(rhs);
    free(obj);
    free(x);
    return (code ? -100 : 0);
}

int LP_NW(EMContext *ctx, int b)
{
    return solve_ballot_simplex(ctx, b, 1);
}
int LPW(EMContext *ctx, int b)
{
    return solve_ballot_simplex(ctx, b, 0);
}
