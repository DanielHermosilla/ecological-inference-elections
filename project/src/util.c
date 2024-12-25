#include "matrixUtils.h"
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// note: Refer to https://www.seehuhn.de/pages/linear.html#blas
// double m[] = {
// 3, 1, 3,
// 1, 5, 9,
// 2, 6, 5
// };

// note: THIS IS THE STRUCT DEFINED IN THE matrixUtils.h file, this comment is just for clarity
// typedef struct
//{
//  double *data; // Pointer to matrix data array (row-major order)
//  int rows;     // Number of rows
//  int cols;     // Number of columns
//} Matrix;

double getInitialP(Matrix x, Matrix w, const char *p_method)
{

    /**
     * @brief Computes the initial probability of the EM algoritm.
     *
     * Given the observables results, it computes a convenient initial "p" value for initiating the
     * algorithm. Currently it supports the "uniform", "group_proportional" and "proportional" methods.
     *
     * @param[in] x Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
     * @param[in] w Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
     * @param[in] p_method The method for calculating the initial parameter. Currently it supports "uniform",
     * "group_proportional" and "proportional" methods.
     * @return Matrix of dimension (gxc) with the initial probability for each demographic group "g" voting for a given
     * candidate "c".
     * @note This should be used only that the first iteration of the EM-algorithm.
     * @warning
     * - Pointers shouldn't be NULL.
     * - `x` and `w` dimensions must be coherent.
     *
     * @see Reference to other relevant functions or documentation.
     */

    // Validation, checks NULL pointer
    if (!x.data || !w.data)
    {
        fprintf(stderr, "A NULL pointer was handed to getInitialP.\n");
        exit(EXIT_FAILURE); // Frees used memory
    }

    // Validation, checks dimentional coherence
    if (x.cols != w.rows && x.cols > 0)
    {
        fprintf(stderr,
                "The dimensions of the matrices handed to getInitialP are incorrect; `x` columns and `w` rows length "
                "must be the same, but they're %d and %d respectivately.\n",
                x.cols, w.rows);
        exit(EXIT_FAILURE);
    }

    const int candidates = x.rows;
    const int groups = w.cols;
    const int ballots = x.cols;

    Matrix *probabilities = (Matrix *)malloc(ballots * sizeof(Matrix)); // Create an allocated array of size ballots.
    free(probabilities); // TODO: When finished, remove this. It's here to remember to free memory.

    // Asumes a Uniform distribution among candidates
    if (strcmp(p_method, "uniform") == 0)
    {
        double baseArray[candidates * groups];

        // Don't parallelize small loops
        for (int i = 0; i < candidates * groups; i++)
        {
            baseArray[i] = 1.0 / candidates;
        }

#pragma omp parallel for
        for (int b = 0; b < ballots; b++)
        {
            // It will fill every .data element with a pointer to baseArray. In this case it won't matter since every
            // element is dependent of each other.
            probabilities[b].rows = candidates;
            probabilities[b].cols = groups;
            probabilities[b].data = baseArray;
        }
    }

    else if (strcmp(p_method, "uniform") == 0)
    {
    }

    double a = 4.0;
    return a;
}

SEXP hello_gsl()
{
    printf("Hello, World from C!\n");

    return R_NilValue;
}
