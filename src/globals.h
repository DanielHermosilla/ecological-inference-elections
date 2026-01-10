// ---- Avoid circular dependencies
#ifndef GLOBALS_H_EIM
#define GLOBALS_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
// Macro for accessing a 3D flattened array (b x g x c)
#define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[(bIdx) * (G) * (C) + (cIdx) * (G) + (gIdx)])
#define MATRIX_AT(matrix, i, j) (matrix.data[(j) * (matrix.rows) + (i)])
#define MATRIX_AT_PTR(matrix, i, j) (matrix->data[(j) * (matrix->rows) + (i)])
#define HSET(ctx, b, g) (&((ctx)->hset[(b) * (ctx)->G + (g)]))
#define KSET(ctx, b, g) (&((ctx)->kset[(b) * (ctx)->G + (g)]))

    // ---- Define the structure to store the input parameters ---- //
    typedef struct
    {
        int S, M;        // Parameters for "MCMC"
        int iters;       // Parameters for importance sampling
        int burnInSteps; // For MCMC and Metropolis
        double stepping_gap;
        char *sampling_method; // For Metropolis
        char *initial_value;   // For Metropolis
        int miniter;
        bool computeLL;
        int monteCarloIter;           // For "MVN CDF"
        double errorThreshold;        // For "MVN CDF"
        const char *simulationMethod; // For "MVN CDF"
        char *prob_cond;
        bool prob_cond_every;
    } QMethodInput;

    // All of the helper functions are made towards double type matrices
    typedef struct
    {
        double *data; // Pointer to matrix data array (col-major order)
        int rows;     // Number of rows
        int cols;     // Number of columns
    } Matrix;

    // The helper functions won't work towards this matrix
    typedef struct
    {
        int *data; // Pointer to matrix data array (col-major order)
        int rows;  // Number of rows
        int cols;  // Number of columns
    } IntMatrix;

    typedef struct
    {
        uint32_t b;
        uint16_t g;
        size_t **data;
        size_t size;
    } Set;

    typedef struct
    {
        uint32_t b;
        IntMatrix *data;
        int *counts;
        size_t size;
    } OmegaSet;

    typedef struct
    {
        // matrices de entrada
        Matrix X; // C×B
        Matrix W; // B×G
        IntMatrix intX;
        IntMatrix intW;
        Matrix probabilities; // G×C
        double *q;            // B×G×C
        double *predicted_votes;
        int iteration;

        Matrix Wnorm;
        double *cdf_seeds;

        // Sizes
        uint32_t B; // TOTAL_BALLOTS
        uint16_t C; // TOTAL_CANDIDATES
        uint16_t G; // TOTAL_GROUPS

        // Precomputation
        uint16_t *ballots_votes;    // length B
        double *inv_ballots_votes;  // length B
        uint32_t *candidates_votes; // length C
        uint32_t *group_votes;      // length G
        double total_votes;
        double *scale_factors; // length B

        // Precomputation
        OmegaSet **omegaset;  // length B
        Set *kset;            // length B
        Set *hset;            // length B
        double **multinomial; // [b][s]
        double *logGamma;     // length max(W)
        double **Qconstant;   // [b][s]

    } EMContext;

    typedef struct
    {
        void (*computeQ)(EMContext *ctx, QMethodInput, double *ll); // Function pointer for computing q
        QMethodInput params;                                        // Holds method-specific parameters
    } QMethodConfig;

    extern uint32_t TOTAL_VOTES;
    extern uint32_t TOTAL_BALLOTS;
    extern uint16_t TOTAL_CANDIDATES;
    extern uint16_t TOTAL_GROUPS;

#ifdef __cplusplus
}
#endif
#endif // GLOBALS_H
