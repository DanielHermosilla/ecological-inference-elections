// GLPK-free dense primal simplex (two-phase), with per-ballot rescale of W if totals mismatch.

#include "globals.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

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
    void *p = Calloc(n * sz, char);
    if (!p)
        error("Allocation error");
    return p;
}
static inline int idxRC(int r, int ld, int c)
{
    return r * ld + c;
}

// ---- Inputs for LP solver ----
typedef struct
{
    const Matrix *X;
    const Matrix *W;
    Matrix *q_bgc;
    double *q;
    int G;
    int C;
    bool x_is_cb;
} LPSolverInput;

static inline double lp_get_x(const LPSolverInput *input, int b, int c)
{
    if (input->x_is_cb)
        return MATRIX_AT_PTR(input->X, c, b);
    return MATRIX_AT_PTR(input->X, b, c);
}
static inline double lp_get_w(const LPSolverInput *input, int b, int g)
{
    return MATRIX_AT_PTR(input->W, b, g);
}
static inline double lp_get_q(const LPSolverInput *input, int b, int g, int c)
{
    if (input->q_bgc != NULL)
        return MATRIX_AT(input->q_bgc[b], g, c);
    return Q_3D(input->q, b, g, c, input->G, input->C);
}
static inline void lp_set_q(const LPSolverInput *input, int b, int g, int c, double value)
{
    if (input->q_bgc != NULL)
        MATRIX_AT(input->q_bgc[b], g, c) = value;
    else
        Q_3D(input->q, b, g, c, input->G, input->C) = value;
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
        Free(P1.T);
        Free(P1.basis);
        return 1;
    }
    SimplexTab P2 = {0};
    build_phase2_from_phase1(&P1, c, &P2);
    Free(P1.T);
    Free(P1.basis);
    rc = simplex_solve_tab(&P2);
    if (rc == 1)
    {
        Free(P2.T);
        Free(P2.basis);
        return 2;
    }
    extract_solution(&P2, x);
    Free(P2.T);
    Free(P2.basis);
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
static LPSolverInput make_lp_input_ctx(EMContext *ctx)
{
    LPSolverInput input = {0};
    input.X = &ctx->X;
    input.W = &ctx->W;
    input.q_bgc = NULL;
    input.q = ctx->q;
    input.G = ctx->G;
    input.C = ctx->C;
    input.x_is_cb = true;
    return input;
}

static LPSolverInput make_lp_input_matrix(const Matrix *X, const Matrix *W, Matrix *q_bgc)
{
    LPSolverInput input = {0};
    input.X = X;
    input.W = W;
    input.q_bgc = q_bgc;
    input.q = NULL;
    input.G = W->cols;
    input.C = X->cols;
    input.x_is_cb = false;
    return input;
}

static int solve_ballot_simplex(const LPSolverInput *input, int b, int mode)
{
    int G = input->G, C = input->C;
    Layout L = mk_layout(G, C);
    double sum_w = 0, sum_x = 0;
    for (int g = 0; g < G; ++g)
        sum_w += lp_get_w(input, b, g);
    for (int c = 0; c < C; ++c)
        sum_x += lp_get_x(input, b, c);
    double rel = (fmax(fabs(sum_w), fabs(sum_x)) > 0) ? fabs(sum_w - sum_x) / fmax(fabs(sum_w), fabs(sum_x)) : 0;
    double alpha = (rel > 1e-12 && sum_w > 0) ? (sum_x / sum_w) : 1.0;
    double *Weff = (double *)xcalloc(G, sizeof(double));
    for (int g = 0; g < G; ++g)
        Weff[g] = lp_get_w(input, b, g) * alpha;

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
            double t = lp_get_q(input, b, g, c);
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
        double xcb = lp_get_x(input, b, c);
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
                lp_set_q(input, b, g, c, 0.0);
    }
    else
    {
        for (int g = 0; g < G; ++g)
            for (int c = 0; c < C; ++c)
            {
                double q = x[col_q(&L, g, c)];
                if (q < 0 && q > -1e-12)
                    q = 0;
                lp_set_q(input, b, g, c, q);
            }
    }
    Free(Weff);
    Free(A);
    Free(rhs);
    Free(obj);
    Free(x);
    return (code ? -100 : 0);
}

int LP_NW(const Matrix *X, const Matrix *W, Matrix *q_bgc, int b)
{
    LPSolverInput input = make_lp_input_matrix(X, W, q_bgc);
    return solve_ballot_simplex(&input, b, 1);
}
int LPW(const Matrix *X, const Matrix *W, Matrix *q_bgc, int b)
{
    LPSolverInput input = make_lp_input_matrix(X, W, q_bgc);
    return solve_ballot_simplex(&input, b, 0);
}

int LP_NW_ctx(EMContext *ctx, int b)
{
    LPSolverInput input = make_lp_input_ctx(ctx);
    return solve_ballot_simplex(&input, b, 1);
}
int LPW_ctx(EMContext *ctx, int b)
{
    LPSolverInput input = make_lp_input_ctx(ctx);
    return solve_ballot_simplex(&input, b, 0);
}
