#include "bootstrap.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

Matrix *iterMat(const Matrix *originalMatrix, const int *indexArr, bool sampleCols, int indexStart)
{
    Matrix *finalMat = malloc(sizeof(Matrix));                            // Allocate space for the Matrix struct
    *finalMat = createMatrix(originalMatrix->rows, originalMatrix->cols); // Allocate data

    for (int j = 0; j < originalMatrix->cols; j++)
    {
        for (int i = 0; i < originalMatrix->rows; i++)
        {
            if (sampleCols)
                MATRIX_AT_PTR(finalMat, i, j) = MATRIX_AT_PTR(originalMatrix, i, indexArr[indexStart + j]);
            else
                MATRIX_AT_PTR(finalMat, i, j) = MATRIX_AT_PTR(originalMatrix, indexArr[indexStart + i], j);
        }
    }
    return finalMat;
}
/**
 *  Returns an array of col-major matrices with bootstrapped matrices.
 *
 * @param[in] xmat The original X array
 * @param[in] wmat The original W array
 * @param[in] bootiter The amount of iterations for bootstrapping
 * @param[in] p_method The method for obtaining the initial probability
 * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
 * and Run", "mult", "mvn_cdf", "mvn_pdf" and "exact" methods.
 * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
 * @param[in] maxIter Integer with a threshold of maximum iterations. Usually it's set to 100.
 * @param[in] maxSeconds Double with the value of the maximum amount of seconds to use.
 * @param[in] verbose Wether to verbose useful outputs.
 * @param[in, out] time The time that the algorithm took.
 * @param[in, out] iterTotal Total amount of iterations.
 * @param[in, out] logLLarr The loglikelihood array
 * @param[in, out] finishing_reason The reason that the algorithm has been stopped. It can either be 0, 1, 2, 3,
 * representing a normal convergence, log likelihood decrease, maximum time reached and maximum iterations reached,
 * respectively.
 *
 *
 * @return An allocated array of size bootiter * TOTAL_BALLOTS that stores matrices.
 */
double *bootstrapA(const Matrix *xmat, const Matrix *wmat, int bootiter, const char *q_method, const char *p_method,
                   const double convergence, const int maxIter, const double maxSeconds, const bool verbose,
                   QMethodInput inputParams)
{

    // ---- Initial variables
    int bdim = wmat->rows;
    int samples = bdim * bootiter;
    int matsize = wmat->cols * xmat->rows;
    // ---- Generate the indices for bootstrap ---- //
    int *indices = Calloc(bdim * bootiter, int);
    GetRNGstate();
    for (int j = 0; j < samples; j++)
    {
        indices[j] = (int)(bdim * unif_rand()) % bdim;
    }
    PutRNGstate();
    // ---...--- //

    // ---- Execute the bootstrap algorithm ---- //
    double *results = Calloc(bootiter * matsize, double);
    for (int i = 0; i < bootiter; i++)
    {
        if (verbose)
            Rprintf("Executing the %dth iteration.\n", i);
        // ---- Declare variables for the current iteration
        Matrix *iterX = iterMat(xmat, indices, true, i * bdim);
        Matrix *iterW = iterMat(wmat, indices, false, i * bdim);
        setParameters(iterX, iterW);
        Matrix iterP = getInitialP(p_method);

        // Declare EM variables, they're not used in this case...
        // It could be useful to yield a mean if the user wants to (logLL mean?)
        double time, logLLarr;
        double *qval = NULL;
        int finishing_reason, iterTotal;
        Matrix resultP = EMAlgoritm(&iterP, q_method, convergence, maxIter, maxSeconds, true, &time, &iterTotal,
                                    &logLLarr, &qval, &finishing_reason, inputParams);
        memcpy(&results[i * matsize], resultP.data, matsize * sizeof(double));

        // ---- Release loop allocated variables ---- //
        // freeMatrix(&iterP);
        printf("first here\n");
        cleanup();
        if (strcmp(q_method, "exact") == 0)
        {
            cleanExact();
        }
        else if (strcmp(q_method, "hnr") == 0)
        {
            cleanHitAndRun();
        }
        printf("second here\n");
        freeMatrix(iterX);
        freeMatrix(iterW);
        printf("third here\n");
        free(iterX);
        free(iterW);
        // ---...--- //
    }
    Free(indices);

    return results;
}
