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

// Macro for making an easier indexation
#define MATRIX_AT(matrix, i, j) (matrix.data[(i) * (matrix.cols) + (j)])

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

    // Validation, checks the method string
    if (strcmp(p_method, "uniform") != 0 && strcmp(p_method, "proportional") != 0 &&
        strcmp(p_method, "group proportional") != 0)
    {
        fprintf(stderr, "The method `%s` passed to getInitialP doesn't exist.\n", p_method);
        exit(EXIT_FAILURE);
    }

    const int candidates = x.rows;
    const int groups = w.cols;
    const int ballots = x.cols;

    Matrix *probabilities =
        (Matrix *)calloc(ballots, sizeof(Matrix)); // Create an allocated array of size ballots with ceros
    free(probabilities); // TODO: When finished, remove this. It's here to remember to free memory.

    // Asumes a Uniform distribution among candidates
    if (strcmp(p_method, "uniform") == 0)
    {
#pragma omp parallel for
        for (int b = 0; b < ballots; b++)
        {
            probabilities[b] = createMatrix(candidates, groups);
            fillMatrix(&probabilities[b], 1.0 / candidates);
        }
    }

    // The proportional method calculates the sum of all the votes per ballot box, and assigns that probability to every
    // demographic group according to ballot box;

    else if (strcmp(p_method, "proportional") == 0)
    {
        double sumBox[ballots];
        double total = 0.0;

        // Step 1: Calculate the total amount of votes per ballot box (n1)
        colSum(&x, sumBox);

        // Step 2: Calculate the total sum (n2)
#pragma omp parallel for // Note that for this case, there shouldn't be a problem of collisions since it's a sum.
        for (int b = 0; b < ballots; b++)
        {
            total += sumBox[b];
        }

        // Step 3: The probabilities for the whole groups and candidates per ballot box is (n1/n2)
#pragma omp parallel for
        for (int b = 0; b < ballots; b++)
        {
            probabilities[b] = createMatrix(candidates, groups);
            fillMatrix(&probabilities[b], sumBox[b] / total);
        }
    }
    else // group proportional
    {
        // Step 1: Calculate each candidate contribution per ballot box

        // Matrix cxb with the contribution of each candidate per ballot
        Matrix candidateContribution = createMatrix(candidates, ballots);
        // Total amount of votes per candidate
        double totalCandidate[ballots];
        rowSum(&x, totalCandidate);

#pragma omp parallel for
        for (int i = 0; i < candidates; i++)
        {
            for (int j = 0; j < ballots; j++)
            {
                if (totalCandidate[i] != 0)
                {
                    MATRIX_AT(totalCandidate, i, j) = MATRIX_AT(totalCandidate, i, j) / totalCandidate[i];
                }
                else
                {
                    MATRIX_AT(totalCandidate, i, j) = 0.0;
                }
            }
        }

        // Step 2: Calculate the total sum (n2)
#pragma omp parallel for // Note that for this case, there shouldn't be a problem of collisions since it's a sum.
        for (int b = 0; b < ballots; b++)
        {
            total += sumBox[b];
        }
    }

    double a = 4.0;
    return a;
}

SEXP hello_gsl()
{
    printf("Hello, World from C!\n");

    return R_NilValue;
}
