#include "globals.h" // tus macros: Q_3D, MATRIX_AT, etc.
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h> // for R_CheckUserInterrupt()
#include <Rinternals.h>
#include <Rmath.h>
#include <glpk.h>
#include <stdlib.h>
#include <string.h>

// Helpers 1-based para índices de columnas (3 bloques de tamaño G*C)
typedef struct
{
    int offX;
    int offQ;
    int offY;
    int G, C;
} ColIdx;

static inline int IDX_X(const ColIdx *L, int g, int c)
{
    return L->offX + (g * L->C + c) + 1; // 1-based
}
static inline int IDX_Q(const ColIdx *L, int g, int c)
{
    return L->offQ + (g * L->C + c) + 1; // 1-based  (q del b fijo)
}
static inline int IDX_Y(const ColIdx *L, int g, int c)
{
    return L->offY + (g * L->C + c) + 1; // 1-based  (y del b fijo)
}

/**
 * Resuelve el LP para la boleta b (0-based) con x_gc, q_bgc, y_bgc.
 * Al finalizar, asigna ctx->q[b,g,c] := x_gc (slice de esa boleta).
 *
 * Devuelve 0 si resuelve sin error de GLPK; !=0 si GLPK falla (pero igual
 * intenta copiar algo coherente si hay valores).
 */
#include <glpk.h>

// Índices 1-based para GLPK
#define IDX_Q(nQ, C, g, c) ((g) * (C) + (c) + 1)
#define IDX_Y(nQ, C, g, c) ((nQ) + (g) * (C) + (c) + 1)

// Acceso matrices columna mayor (como usas en tu proyecto)
#define MAT_AT(M, row, col) MATRIX_AT(M, row, col) // alias por claridad
//
int LP_NW(EMContext *ctx, int b)
{
    const int B = (int)ctx->B;
    const int G = (int)ctx->G;
    const int C = (int)ctx->C;

    if (b < 0 || b >= B)
    {
        Rprintf("fastei_glpk_fit_q_single: indice b fuera de rango.\n");
        return -1;
    }

    // Columnas: primero todos q_{gc}, luego todos y_{gc}
    const int nQ = G * C;
    const int nY = G * C;
    const int nCols = nQ + nY;

    // Filas:
    //  - ABS: 2 por cada (g,c)
    //  - Grupo: G filas (sum_c q_{gc} = 1)
    //  - Candidato: C filas (sum_g w_{bg} q_{gc} = x_{cb})
    const int nAbs = 2 * G * C;
    const int nGrp = G;
    const int nCat = C;
    const int nRows = nAbs + nGrp + nCat;

    // No ceros:
    //  - ABS: cada fila involucra (y,q) ⇒ 2 nnz por fila ⇒ 2*(2*G*C) = 4*G*C
    //  - Grupo: cada fila tiene C coef de q ⇒ G*C
    //  - Candidato: cada fila tiene G coef de q ⇒ C*G
    const int nnz = 4 * G * C + G * C + C * G; // = 6*G*C

    glp_prob *lp = glp_create_prob();
    glp_set_prob_name(lp, "LP_qy_one_ballot");
    glp_set_obj_dir(lp, GLP_MIN);

    glp_add_rows(lp, nRows);
    glp_add_cols(lp, nCols);

    // --- Definir columnas: q_{gc} >= 0 (no contribuye al objetivo) ---
    for (int g = 0; g < G; ++g)
    {
        for (int c = 0; c < C; ++c)
        {
            int j = IDX_Q(nQ, C, g, c);
            char name[48];
            snprintf(name, sizeof(name), "q_%d_%d_b%d", g, c, b);
            glp_set_col_name(lp, j, name);
            glp_set_col_bnds(lp, j, GLP_LO, 0.0, 0.0);
            glp_set_obj_coef(lp, j, 0.0);
        }
    }

    // --- Definir columnas: y_{gc} >= 0 (sí entra al objetivo) ---
    for (int g = 0; g < G; ++g)
    {
        for (int c = 0; c < C; ++c)
        {
            int j = IDX_Y(nQ, C, g, c);
            char name[48];
            snprintf(name, sizeof(name), "y_%d_%d_b%d", g, c, b);
            glp_set_col_name(lp, j, name);
            glp_set_col_bnds(lp, j, GLP_LO, 0.0, 0.0);
            glp_set_obj_coef(lp, j, 1.0);
        }
    }

    int *ia = (int *)malloc((nnz + 1) * sizeof(int));
    int *ja = (int *)malloc((nnz + 1) * sizeof(int));
    double *ar = (double *)malloc((nnz + 1) * sizeof(double));
    if (!ia || !ja || !ar)
    {
        if (ia)
            free(ia);
        if (ja)
            free(ja);
        if (ar)
            free(ar);
        glp_delete_prob(lp);
        Rprintf("LP: sin memoria para matriz dispersa.\n");
        return -2;
    }

    int r = 1; // fila actual (1-based)
    int k = 1; // puntero en (ia,ja,ar) 1-based

    // --- (1) Restricciones ABS (2 por cada g,c) ---
    // y_{gc} - w_{bg} q_{gc} >= - w_{bg} t_{gc}
    // y_{gc} + w_{bg} q_{gc} >=   w_{bg} t_{gc}
    for (int g = 0; g < G; ++g)
    {
        const double wbg = MAT_AT(ctx->W, b, g); // W es BxG
        for (int c = 0; c < C; ++c)
        {
            const double tgc =
                Q_3D(ctx->q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES); // MAT_AT(ctx->q, g, c); // parámetro t_{gc}
            const int jq = IDX_Q(nQ, C, g, c);
            const int jy = IDX_Y(nQ, C, g, c);

            // abs1
            {
                char name[48];
                snprintf(name, sizeof(name), "abs1_g%d_c%d_b%d", g, c, b);
                glp_set_row_name(lp, r, name);
                glp_set_row_bnds(lp, r, GLP_LO, -wbg * tgc, 0.0);

                ia[k] = r;
                ja[k] = jy;
                ar[k] = 1.0;
                ++k; // + y
                ia[k] = r;
                ja[k] = jq;
                ar[k] = -wbg;
                ++k; // - w*q
                ++r;
            }

            // abs2
            {
                char name[48];
                snprintf(name, sizeof(name), "abs2_g%d_c%d_b%d", g, c, b);
                glp_set_row_name(lp, r, name);
                glp_set_row_bnds(lp, r, GLP_LO, wbg * tgc, 0.0);

                ia[k] = r;
                ja[k] = jy;
                ar[k] = 1.0;
                ++k; // + y
                ia[k] = r;
                ja[k] = jq;
                ar[k] = wbg;
                ++k; // + w*q
                ++r;
            }
        }
    }

    // --- (2) Suma por grupo: sum_c q_{gc} = 1  (para cada g) ---
    for (int g = 0; g < G; ++g)
    {
        char name[48];
        snprintf(name, sizeof(name), "sum_group_g%d_b%d", g, b);
        glp_set_row_name(lp, r, name);
        glp_set_row_bnds(lp, r, GLP_FX, 1.0, 1.0);

        for (int c = 0; c < C; ++c)
        {
            int jq = IDX_Q(nQ, C, g, c);
            ia[k] = r;
            ja[k] = jq;
            ar[k] = 1.0;
            ++k;
        }
        ++r;
    }

    // --- (3) Suma por candidato: sum_g w_{bg} q_{gc} = x_{cb}  (para cada c) ---
    for (int c = 0; c < C; ++c)
    {
        const double xcb = MAT_AT(ctx->X, c, b); // X es CxB (c,b)
        char name[48];
        snprintf(name, sizeof(name), "sum_cand_c%d_b%d", c, b);
        glp_set_row_name(lp, r, name);
        glp_set_row_bnds(lp, r, GLP_FX, xcb, xcb);

        for (int g = 0; g < G; ++g)
        {
            const double wbg = MAT_AT(ctx->W, b, g);
            int jq = IDX_Q(nQ, C, g, c);
            ia[k] = r;
            ja[k] = jq;
            ar[k] = wbg;
            ++k;
        }
        ++r;
    }

    // Cargar matriz
    glp_load_matrix(lp, nnz, ia, ja, ar);

    // Resolver
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_ERR; // silencioso
    int glp_ret = glp_simplex(lp, &parm);

    // Extraer solución q_{gc} -> ctx->q[b,g,c]
    if (glp_ret == 0)
    {
        // Rprintf("funcionó el LP!\n");
        for (int g = 0; g < G; ++g)
        {
            for (int c = 0; c < C; ++c)
            {
                int jq = IDX_Q(nQ, C, g, c);
                double qgc = glp_get_col_prim(lp, jq);
                if (qgc < 0.0)
                    qgc = 0.0; // limpieza numérica
                Q_3D(ctx->q, b, g, c, G, C) = qgc;
            }
        }
    }
    else
    {
        // fallback: deja q de esta boleta en 0 si falla
        for (int g = 0; g < G; ++g)
            for (int c = 0; c < C; ++c)
                Q_3D(ctx->q, b, g, c, G, C) = 0.0;
    }

    free(ia);
    free(ja);
    free(ar);
    glp_delete_prob(lp);
    return glp_ret; // 0 = OK
}

// Índices 1-based para GLPK (primero todos q, luego todos y)

// Alias para acceso a matrices (col-major) como en tu proyecto

int LPW(EMContext *ctx, int b)
{
    const int B = (int)ctx->B;
    const int G = (int)ctx->G;
    const int C = (int)ctx->C;

    if (b < 0 || b >= B)
    {
        Rprintf("fastei_glpk_fit_q_single: indice b fuera de rango.\n");
        return -1;
    }

    // Columnas: q_{gc} (G*C) y y_{gc} (G*C)
    const int nQ = G * C;
    const int nY = G * C;
    const int nCols = nQ + nY;

    // Filas:
    //  - ABS: 2 por cada (g,c)  => 2*G*C
    //  - Grupo: sum_c q_{gc} = 1   (G filas)
    //  - Candidato: sum_g w_{bg} q_{gc} = x_{cb}  (C filas)
    const int nAbs = 2 * G * C;
    const int nGrp = G;
    const int nCat = C;
    const int nRows = nAbs + nGrp + nCat;

    // No ceros:
    //  - ABS: cada fila involucra y y q -> 2 nnz por fila => 2*(2*G*C) = 4*G*C
    //  - Grupo: G filas * C nnz = G*C
    //  - Candidato: C filas * G nnz = C*G
    const int nnz = 4 * G * C + G * C + C * G; // 6*G*C

    glp_prob *lp = glp_create_prob();
    glp_set_prob_name(lp, "LP_qy_one_ballot");
    glp_set_obj_dir(lp, GLP_MIN);

    glp_add_rows(lp, nRows);
    glp_add_cols(lp, nCols);

    // --- Definir columnas: q_{gc} >= 0 (no entra al objetivo) ---
    for (int g = 0; g < G; ++g)
    {
        for (int c = 0; c < C; ++c)
        {
            int j = IDX_Q(nQ, C, g, c);
            char name[48];
            snprintf(name, sizeof(name), "q_%d_%d_b%d", g, c, b);
            glp_set_col_name(lp, j, name);
            glp_set_col_bnds(lp, j, GLP_LO, 0.0, 0.0);
            glp_set_obj_coef(lp, j, 0.0);
        }
    }

    // --- Definir columnas: y_{gc} >= 0 (sí entra al objetivo) ---
    for (int g = 0; g < G; ++g)
    {
        for (int c = 0; c < C; ++c)
        {
            int j = IDX_Y(nQ, C, g, c);
            char name[48];
            snprintf(name, sizeof(name), "y_%d_%d_b%d", g, c, b);
            glp_set_col_name(lp, j, name);
            glp_set_col_bnds(lp, j, GLP_LO, 0.0, 0.0);
            glp_set_obj_coef(lp, j, 1.0); // min sum y
        }
    }

    int *ia = (int *)malloc((nnz + 1) * sizeof(int));
    int *ja = (int *)malloc((nnz + 1) * sizeof(int));
    double *ar = (double *)malloc((nnz + 1) * sizeof(double));
    if (!ia || !ja || !ar)
    {
        if (ia)
            free(ia);
        if (ja)
            free(ja);
        if (ar)
            free(ar);
        glp_delete_prob(lp);
        Rprintf("LP: sin memoria para matriz dispersa.\n");
        return -2;
    }

    int r = 1; // fila actual (1-based)
    int k = 1; // puntero en (ia,ja,ar) 1-based

    // --- (1) Restricciones ABS (2 por cada g,c) ---
    // y_{gc} >=  q_{gc} - t_{gc}  <=>  y - q >= -t
    // y_{gc} >= -q_{gc} + t_{gc}  <=>  y + q >=  t
    for (int g = 0; g < G; ++g)
    {
        for (int c = 0; c < C; ++c)
        {

            const double tgc =
                Q_3D(ctx->q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES); // MAT_AT(ctx->q, g, c); // parámetro t_{gc}
            const int jq = IDX_Q(nQ, C, g, c);
            const int jy = IDX_Y(nQ, C, g, c);

            // y - q >= -t
            {
                char name[48];
                snprintf(name, sizeof(name), "abs1_g%d_c%d_b%d", g, c, b);
                glp_set_row_name(lp, r, name);
                glp_set_row_bnds(lp, r, GLP_LO, -tgc, 0.0);

                ia[k] = r;
                ja[k] = jy;
                ar[k] = 1.0;
                ++k; // + y
                ia[k] = r;
                ja[k] = jq;
                ar[k] = -1.0;
                ++k; // - q
                ++r;
            }

            // y + q >=  t
            {
                char name[48];
                snprintf(name, sizeof(name), "abs2_g%d_c%d_b%d", g, c, b);
                glp_set_row_name(lp, r, name);
                glp_set_row_bnds(lp, r, GLP_LO, tgc, 0.0);

                ia[k] = r;
                ja[k] = jy;
                ar[k] = 1.0;
                ++k; // + y
                ia[k] = r;
                ja[k] = jq;
                ar[k] = 1.0;
                ++k; // + q
                ++r;
            }
        }
    }

    // --- (2) Suma por grupo: sum_c q_{gc} = 1 (∀g) ---
    for (int g = 0; g < G; ++g)
    {
        char name[48];
        snprintf(name, sizeof(name), "sum_group_g%d_b%d", g, b);
        glp_set_row_name(lp, r, name);
        glp_set_row_bnds(lp, r, GLP_FX, 1.0, 1.0);

        for (int c = 0; c < C; ++c)
        {
            int jq = IDX_Q(nQ, C, g, c);
            ia[k] = r;
            ja[k] = jq;
            ar[k] = 1.0;
            ++k;
        }
        ++r;
    }

    // --- (3) Suma por candidato: sum_g w_{bg} q_{gc} = x_{cb} (∀c) ---
    for (int c = 0; c < C; ++c)
    {
        const double xcb = MAT_AT(ctx->X, c, b); // X es CxB (c,b)
        char name[48];
        snprintf(name, sizeof(name), "sum_cand_c%d_b%d", c, b);
        glp_set_row_name(lp, r, name);
        glp_set_row_bnds(lp, r, GLP_FX, xcb, xcb);

        for (int g = 0; g < G; ++g)
        {
            const double wbg = MAT_AT(ctx->W, b, g);
            int jq = IDX_Q(nQ, C, g, c);
            ia[k] = r;
            ja[k] = jq;
            ar[k] = wbg;
            ++k;
        }
        ++r;
    }

    // Cargar matriz
    glp_load_matrix(lp, nnz, ia, ja, ar);

    // Resolver
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_ERR; // silencioso (solo errores)
    int glp_ret = glp_simplex(lp, &parm);

    // Extraer solución: q_{gc} -> ctx->q[b,g,c]
    if (glp_ret == 0)
    {
        for (int g = 0; g < G; ++g)
        {
            for (int c = 0; c < C; ++c)
            {
                int jq = IDX_Q(nQ, C, g, c);
                double qgc = glp_get_col_prim(lp, jq);
                if (qgc < 0.0)
                    qgc = 0.0; // limpieza numérica
                Q_3D(ctx->q, b, g, c, G, C) = qgc;
            }
        }
    }
    else
    {
        // fallback: deja q de esta boleta en 0 si falla
        for (int g = 0; g < G; ++g)
            for (int c = 0; c < C; ++c)
                Q_3D(ctx->q, b, g, c, G, C) = 0.0;
    }

    free(ia);
    free(ja);
    free(ar);
    glp_delete_prob(lp);
    return glp_ret; // 0 = OK
}
