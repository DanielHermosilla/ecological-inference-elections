// ---- Avoid circular dependencies
#ifndef GLOBALS_H_EIM
#define GLOBALS_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include <stdint.h>
#include <stdio.h>

// Macro for accessing a 3D flattened array (b x g x c)
// #define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[((bIdx) * (G) * (C)) + ((gIdx) * (C)) + (cIdx)])
#define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[(bIdx) * (G) * (C) + (cIdx) * (G) + (gIdx)])
#define MATRIX_AT(matrix, i, j) (matrix.data[(j) * (matrix.rows) + (i)])
#define MATRIX_AT_PTR(matrix, i, j) (matrix->data[(j) * (matrix->rows) + (i)])

    // ---- Define the structure to store the input parameters ---- //
    typedef struct
    {
        int S, M;                     // Parameters for ""Hit and Run"
        int monteCarloIter;           // For "MVN CDF"
        double errorThreshold;        // For "MVN CDF"
        const char *simulationMethod; // For "MVN CDF"
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
        size_t *data; // Pointer to matrix data array (col-major order)
        int rows;     // Number of rows
        int cols;     // Number of columns
    } SizeTMatrix;

    extern uint32_t TOTAL_VOTES;
    extern uint32_t TOTAL_BALLOTS;
    extern uint16_t TOTAL_CANDIDATES;
    extern uint16_t TOTAL_GROUPS;
    extern uint16_t *BALLOTS_VOTES;    // Total votes per ballot
    extern uint32_t *CANDIDATES_VOTES; // Total votes per candidate
    extern uint32_t *GROUP_VOTES;      // Total votes per group
    extern double *inv_BALLOTS_VOTES;
    extern Matrix *X;
    extern Matrix *W;
#ifdef __cplusplus
}
#endif
#endif // GLOBALS_H
