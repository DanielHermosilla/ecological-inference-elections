// GLPK-free, dense two-phase primal simplex using BLAS.
//
// Vars (all >= 0):
//   q_{g,c}  (G*C)
//   y_{g,c}  (G*C)
// Internals we add for standard form (not in your objective):
//   s1_{g,c}, s2_{g,c}  (slacks for the two ABS ≥-constraints)
//
// Objective (min):
//   sum_{g,c} y_{g,c}
//
// Constraints:
//   ABS (weighted or unweighted via mode):
//     y_{g,c} - coef*q_{g,c} >= -coef*t_{g,c}
//     y_{g,c} + coef*q_{g,c} >=  coef*t_{g,c}
//   GROUP (∀g):
//     sum_c q_{g,c} = 1
//   CAND  (∀c):
//     sum_g w_{b,g} q_{g,c} = X_{c,b}
//
// Public API preserved:
//   int LP_NW(EMContext *ctx, int b);  // weighted ABS (same as your GLPK LP_NW)
//   int LPW  (EMContext *ctx, int b);  // unweighted ABS (same as your GLPK LPW)

#include "globals.h" // EMContext, MATRIX_AT, Q_3D, TOTAL_* macros
#include <R_ext/BLAS.h>
#include <R_ext/RS.h> // F77_CALL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// -------- BLAS row helpers --------
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

// -------- utils --------
static void *xcalloc(size_t n, size_t sz)
{
    void *p = calloc(n, sz);
    if (!p)
    {
        perror("calloc");
        exit(1);
    }
    return p;
}
static inline int idxRC(int r, int ld, int c)
{
    return r * ld + c;
} // row-major, ld = nCols

// -------- simplex tableau --------
typedef struct
{
    int m;      // #rows (constraints)
    int n;      // #cols (vars, excl RHS)
    int ncols;  // n + 1 (RHS)
    double *T;  // (m+1) x (n+1) tableau; last row = objective; last col = RHS
    int *basis; // size m; basis column index (0..n-1) or -1
    double tol;
} SimplexTab;

static void pivot_at(SimplexTab *tab, int leave, int enter)
{
    const int m = tab->m, ncols = tab->ncols;
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

// Minimization: enter any column with NEGATIVE reduced cost (Bland)
static int choose_entering_bland(const SimplexTab *tab)
{
    const double *obj = &tab->T[idxRC(tab->m, tab->ncols, 0)];
    for (int j = 0; j < tab->n; ++j)
        if (obj[j] < -tab->tol)
            return j;
    return -1; // optimal
}

static int choose_leaving(const SimplexTab *tab, int enter)
{
    const int m = tab->m, ncols = tab->ncols;
    const double *T = tab->T;
    int leave = -1;
    double best = HUGE_VAL;
    for (int i = 0; i < m; ++i)
    {
        double a = T[idxRC(i, ncols, enter)];
        if (a > tab->tol)
        {
            double rhs = T[idxRC(i, ncols, tab->n)];
            double ratio = rhs / a;
            if (ratio < best - 1e-12 ||
                (fabs(ratio - best) <= 1e-12 && tab->basis[i] > (leave >= 0 ? tab->basis[leave] : -1)))
            { // Bland-ish tie-break
                best = ratio;
                leave = i;
            }
        }
    }
    return leave; // -1 => unbounded
}

static int simplex_solve_tab(SimplexTab *tab)
{
    while (1)
    {
        int e = choose_entering_bland(tab);
        if (e < 0)
            return 0; // optimal
        int l = choose_leaving(tab, e);
        if (l < 0)
            return 1; // unbounded
        pivot_at(tab, l, e);
    }
}

// Phase I builder: A(mxn), b(m)  -> add m artificials, initial basis = artificials
static void build_phase1_full(const double *A, const double *b, int m, int n, SimplexTab *tab)
{
    const int n_art = m;
    tab->m = m;
    tab->n = n + n_art;
    tab->ncols = tab->n + 1;
    tab->T = (double *)xcalloc((size_t)(m + 1) * (size_t)(tab->ncols), sizeof(double));
    tab->basis = (int *)xcalloc(m, sizeof(int));
    tab->tol = 1e-10;

    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
            tab->T[idxRC(i, tab->ncols, j)] = A[idxRC(i, n, j)];
        int ja = n + i;
        tab->T[idxRC(i, tab->ncols, ja)] = 1.0;      // artificial
        tab->T[idxRC(i, tab->ncols, tab->n)] = b[i]; // RHS
        tab->basis[i] = ja;
    }

    // Objective: minimize sum artificials => obj = 1 on artificials
    double *obj = &tab->T[idxRC(m, tab->ncols, 0)];
    for (int j = n; j < n + n_art; ++j)
        obj[j] = 1.0;

    // Canonicalize (reduced costs) wrt current basic artificials
    for (int i = 0; i < m; ++i)
    {
        double *row = &tab->T[idxRC(i, tab->ncols, 0)];
        row_axpy(tab->ncols, -1.0, row, obj); // make artificial reduced costs 0
    }
}

// Sum of artificial basic values (always >=0 after row-normalization). If > eps -> infeasible.
static double phase1_artificial_value(const SimplexTab *p1)
{
    const int n_orig = p1->n - p1->m;
    const int ncols = p1->ncols;
    double s = 0.0;
    for (int i = 0; i < p1->m; ++i)
    {
        int b = p1->basis[i];
        if (b >= n_orig)
            s += p1->T[idxRC(i, ncols, p1->n)];
    }
    return s;
}

// Pivot-out basic artificials if there is support in original columns.
static int phase1_eliminate_artificials(SimplexTab *p1)
{
    const int m = p1->m;
    const int n_orig = (p1->n) - m;
    const int ncols = p1->ncols;

    for (int i = 0; i < m; ++i)
    {
        int b = p1->basis[i];
        if (b >= n_orig)
        { // artificial basic
            int enter = -1;
            double best = 0.0;
            for (int j = 0; j < n_orig; ++j)
            {
                double aij = p1->T[idxRC(i, ncols, j)];
                if (fabs(aij) > best + 1e-12)
                {
                    best = fabs(aij);
                    enter = j;
                }
            }
            if (enter >= 0)
            {
                pivot_at(p1, i, enter);
            }
            else
            {
                // No original support; row must be redundant, ensure RHS ~ 0
                double rhs = p1->T[idxRC(i, ncols, p1->n)];
                if (fabs(rhs) > 1e-9)
                    return 1; // inconsistent -> infeasible
            }
        }
    }
    return 0;
}

static void set_phase2_objective(SimplexTab *tab, const double *c, int n_orig)
{
    double *obj = &tab->T[idxRC(tab->m, tab->ncols, 0)];
    for (int j = 0; j < tab->n; ++j)
        obj[j] = (j < n_orig ? c[j] : 0.0);
    obj[tab->n] = 0.0;
    // canonicalize vs current basis
    for (int i = 0; i < tab->m; ++i)
    {
        int b = tab->basis[i];
        if (b >= 0)
        {
            double coef = obj[b];
            if (fabs(coef) > 0)
            {
                double *row = &tab->T[idxRC(i, tab->ncols, 0)];
                row_axpy(tab->ncols, -coef, row, obj);
            }
        }
    }
}

// Compact P1 -> P2 keep only original columns
static void build_phase2_from_phase1(SimplexTab *p1, const double *c, SimplexTab *p2)
{
    const int m = p1->m;
    const int n_orig = (p1->n) - m;

    // 1) try to remove artificials from basis; if impossible and RHS ≠ 0 => infeasible
    if (phase1_eliminate_artificials(p1))
    {
        // mark infeasible
        p1->T[idxRC(m, p1->ncols, p1->n)] = HUGE_VAL;
    }

    // 2) set Phase-II objective on p1, canonicalize
    set_phase2_objective(p1, c, n_orig);

    // 3) compact to P2 (original columns only)
    p2->m = m;
    p2->n = n_orig;
    p2->ncols = n_orig + 1;
    p2->T = (double *)xcalloc((size_t)(m + 1) * (size_t)(p2->ncols), sizeof(double));
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
    p2->T[idxRC(m, p2->ncols, p2->n)] = p1->T[idxRC(m, p1->ncols, p1->n)];

    // Re-establish clean basis on P2 (pivot row-normal if needed)
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

static void extract_solution(const SimplexTab *p2, double *x /*size n*/)
{
    const int m = p2->m, n = p2->n, ncols = p2->ncols;
    memset(x, 0, sizeof(double) * n);
    for (int i = 0; i < m; ++i)
    {
        int b = p2->basis[i];
        if (0 <= b && b < n)
            x[b] = p2->T[idxRC(i, ncols, n)];
    }
}

// Solve standard-form equalities with x >= 0
// returns: 0 ok, 1 infeasible, 2 unbounded
static int simplex_solve_dense_equalities(double *A, double *b, const double *c, int m, int n, double *x_out)
{
    // Pre: ensure b >= 0 so Phase I initial BFS (artificial = b) is feasible
    for (int i = 0; i < m; ++i)
    {
        if (b[i] < 0)
        {
            row_scale(n, &A[idxRC(i, n, 0)], -1.0);
            b[i] = -b[i];
        }
    }

    // Phase I
    SimplexTab P1 = {0};
    build_phase1_full(A, b, m, n, &P1);
    int rc = simplex_solve_tab(&P1);
    double art_val = phase1_artificial_value(&P1);
    if (rc == 1 || art_val > 1e-7)
    {
        free(P1.T);
        free(P1.basis);
        return 1; // infeasible
    }

    // Phase II
    SimplexTab P2 = {0};
    build_phase2_from_phase1(&P1, c, &P2);
    double rhs_obj = P2.T[idxRC(P2.m, P2.ncols, P2.n)];
    free(P1.T);
    free(P1.basis);
    if (!isfinite(rhs_obj))
    {
        free(P2.T);
        free(P2.basis);
        return 1;
    }

    rc = simplex_solve_tab(&P2);
    if (rc == 1)
    {
        free(P2.T);
        free(P2.basis);
        return 2;
    } // unbounded

    extract_solution(&P2, x_out);
    free(P2.T);
    free(P2.basis);
    return 0;
}

// -------- Build (q,y,s) model --------
typedef struct
{
    int G, C;
    int nQ, nY, nS, nVars;
    int nAbs, nGrp, nCat, nRows;
    int offQ, offY, offS;
} Layout;

static Layout mk_layout(int G, int C)
{
    Layout L;
    L.G = G;
    L.C = C;
    L.nQ = G * C;
    L.nY = G * C;
    L.nAbs = 2 * G * C;
    L.nS = L.nAbs; // one slack per ABS inequality row
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

// mode == 1  -> weighted ABS (GLPK LP_NW)
// mode == 0  -> unweighted ABS (GLPK LPW)
static int solve_ballot_simplex(EMContext *ctx, int ballot, int mode)
{
    const int G = (int)ctx->G, C = (int)ctx->C;
    Layout L = mk_layout(G, C);

    // Build dense A (m x n), b (m), c (n)
    double *A = (double *)xcalloc((size_t)L.nRows * (size_t)L.nVars, sizeof(double));
    double *rhs = (double *)xcalloc(L.nRows, sizeof(double));
    double *obj = (double *)xcalloc(L.nVars, sizeof(double));

    // Objective: min sum y
    for (int g = 0; g < G; ++g)
        for (int c = 0; c < C; ++c)
            obj[col_y(&L, g, c)] = 1.0;

    // ABS rows (exactly as your GLPK, with slacks to make equalities)
    for (int g = 0; g < G; ++g)
    {
        const double wbg = MATRIX_AT(ctx->W, ballot, g); // W is BxG
        for (int c = 0; c < C; ++c)
        {
            const double tgc = Q_3D(ctx->q, ballot, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES);
            const double coef = (mode ? wbg : 1.0);

            // y - coef*q - s1 = -coef*t
            {
                int r = row_abs1(&L, g, c);
                A[idxRC(r, L.nVars, col_y(&L, g, c))] = 1.0;
                A[idxRC(r, L.nVars, col_q(&L, g, c))] = -coef;
                A[idxRC(r, L.nVars, col_s1(&L, g, c))] = -1.0;
                rhs[r] = -coef * tgc;
            }
            // y + coef*q - s2 =  coef*t
            {
                int r = row_abs2(&L, g, c);
                A[idxRC(r, L.nVars, col_y(&L, g, c))] = 1.0;
                A[idxRC(r, L.nVars, col_q(&L, g, c))] = coef;
                A[idxRC(r, L.nVars, col_s2(&L, g, c))] = -1.0;
                rhs[r] = coef * tgc;
            }
        }
    }

    // GROUP rows: sum_c q_{g,c} = 1
    for (int g = 0; g < G; ++g)
    {
        int r = row_grp(&L, g);
        for (int c = 0; c < C; ++c)
            A[idxRC(r, L.nVars, col_q(&L, g, c))] = 1.0;
        rhs[r] = 1.0;
    }

    // CAND rows: sum_g w_{b,g} q_{g,c} = X_{c,b}
    for (int c = 0; c < C; ++c)
    {
        int r = row_cand(&L, c);
        double xcb = MATRIX_AT(ctx->X, c, ballot); // X is CxB
        for (int g = 0; g < G; ++g)
        {
            double wbg = MATRIX_AT(ctx->W, ballot, g);
            A[idxRC(r, L.nVars, col_q(&L, g, c))] = wbg;
        }
        rhs[r] = xcb;
    }

    // Solve in standard form
    double *x = (double *)xcalloc(L.nVars, sizeof(double));
    int code = simplex_solve_dense_equalities(A, rhs, obj, L.nRows, L.nVars, x);

    int rc = 0;
    if (code != 0)
    {
        // GLPK fallback behavior: zero out this ballot’s q if infeasible/unbounded
        for (int g = 0; g < G; ++g)
            for (int c = 0; c < C; ++c)
                Q_3D(ctx->q, ballot, g, c, G, C) = 0.0;
        rc = (code == 1 ? -100 : -101);
    }
    else
    {
        // copy q back (and clamp tiny negatives like GLPK path)
        for (int g = 0; g < G; ++g)
            for (int c = 0; c < C; ++c)
            {
                double qgc = x[col_q(&L, g, c)];
                if (qgc < 0.0)
                    qgc = 0.0; // numeric cleanup (GLPK also clamps)
                Q_3D(ctx->q, ballot, g, c, G, C) = qgc;
            }
    }

    free(A);
    free(rhs);
    free(obj);
    free(x);
    return rc;
}

// ---- Public wrappers (names & semantics match your GLPK code) ----
int LP_NW(EMContext *ctx, int b)
{ // weighted ABS (GLPK LP_NW)
    return solve_ballot_simplex(ctx, b, /*mode=*/1);
}
int LPW(EMContext *ctx, int b)
{ // unweighted ABS (GLPK LPW)
    return solve_ballot_simplex(ctx, b, /*mode=*/0);
}
