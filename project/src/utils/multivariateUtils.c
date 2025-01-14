#include "multivariateUtils.h"
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Computes the parameters of the unconditional probability without the last candidate
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution.
 *
 * @param[in] b The index of the ballot box
 * @param[in] g The index of the group
 * @param[in] *probabilitiesReduce Matrix of dimension (gxc-1) with the probabilities of each group and candidate,
 * except the last one. Consider that the information of the last candidate is redundant. The probability matrix could
 * be reduced with the "removeRows()" function from matrixUtils.
 * @param[in, out] *mu An array of size c-1 to store the results of the average.
 * @param[in, out] *sigma A matrix of size (c-1, c-1) to store the sigma matrix.
 *
 * @warning Remember to eliminate one dimension for candidates.
 *
 * @return void. Results to be written on mu and sigma.
 *
 */

void getParams(int b, const Matrix *probabilitiesReduced, double *mu, Matrix *sigma)
{

    // ---- Check parameters ---- //
    // ---- Note: Mu must be of size TOTAL_CANDIDATES-1 and sigma of size (TOTAL_CANDIDATES-1xTOTAL_CANDIDATES-1) ----
    if (probabilitiesReduced->cols != TOTAL_CANDIDATES - 1)
    {
        fprintf(
            stderr,
            "The probability matrix handed should consider C-1 candidates, but it has %d columns. Consider using the "
            "`removeLastCols()` function from matrixUtils.\n",
            probabilitiesReduced->cols);
        exit(EXIT_FAILURE);
    }
    if (W == NULL && X == NULL)
    {
        fprintf(stderr, "The `w` and `x` matrices aren't defined.\n");
        exit(EXIT_FAILURE);
    }
    // --- ... --- //

    // --- Calculations --- //
    // ---- The votes that a group has made on a given ballot ----
    double *groupVotesPerBallot = getRow(W, b);

    // ---- Computation of mu ----
    // ---- Performing the matrix multiplication of p^T * w_b
    cblas_dgemv(CblasRowMajor,              // Row-major order
                CblasTrans,                 // Transpose the matrix (p^T)
                TOTAL_GROUPS,               // Rows of original matrix (G)
                TOTAL_CANDIDATES - 1,       // Columns of original matrix (C)
                1.0,                        // alpha
                probabilitiesReduced->data, // Matrix p
                TOTAL_CANDIDATES - 1,       // Leading dimension (rows in original p)
                groupVotesPerBallot,        // Vector w_b
                1,                          // Increment for x
                0.0,                        // beta
                mu,                         // Output vector μ
                1                           // Increment for y
    );

    // ---- Computation of sigma ----
    // ---- Get a diagonal matrix with the group votes on a given ballot ----
    Matrix diagonalVotesPerBallot = createDiagonalMatrix(groupVotesPerBallot, TOTAL_GROUPS);
    // ---- Temporary matrix to store results ----
    Matrix temp = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_GROUPS);
    // ---- Calculates the matrix multiplication of p^T * diag(w_b); result must be (C-1 x G) ----
    cblas_dgemm(CblasRowMajor,               // Row-major order
                CblasTrans,                  // No transpose (Matrix A => p^T)
                CblasNoTrans,                // No transpose
                TOTAL_CANDIDATES - 1,        // Rows
                TOTAL_GROUPS,                // Columns
                TOTAL_GROUPS,                // Inner dimension
                1.0,                         // alpha
                probabilitiesReduced->data,  // Matrix p
                TOTAL_CANDIDATES - 1,        // Leading dimension
                diagonalVotesPerBallot.data, // Matrix diag(w_b)
                TOTAL_GROUPS,                // Leading dimension
                0.0,                         // beta
                temp.data,                   // Output matrix
                TOTAL_GROUPS);               // Leading dimension

    // ---- Calculates the matrix multiplication of (p^T * diag(w_b)) * p; result must be (C-1 x C-1) ----
    cblas_dgemm(CblasRowMajor,              // Row-major order
                CblasNoTrans,               // No transpose
                CblasNoTrans,               // No transpose
                TOTAL_CANDIDATES - 1,       // Rows
                TOTAL_CANDIDATES - 1,       // Columns
                TOTAL_GROUPS,               // Inner dimension
                1.0,                        // alpha
                temp.data,                  // Matrix (p^T * diag(w_b))
                TOTAL_GROUPS,               // Leading dimension
                probabilitiesReduced->data, // Matrix p
                TOTAL_CANDIDATES - 1,       // Leading dimension
                0.0,                        // beta
                sigma->data,                // Output Σ_b
                TOTAL_CANDIDATES - 1);      // Leading dimension

    // ---- Substract the diagonal with the average ----
    // ---- Note: This could be optimized with a cBLAS call too ----

    for (int j = 0; j < TOTAL_CANDIDATES - 1; j++)
    { // ---- For each candidate
        for (int i = 0; i < TOTAL_CANDIDATES - 1; i++)
        { // ---- For each candidate given another candidate
            if (i == j)
            { // ---- If it corresponds to a diagonal, substract diagonal
                MATRIX_AT_PTR(sigma, i, j) = mu[i] - MATRIX_AT_PTR(sigma, i, j);
                continue;
            }
            MATRIX_AT_PTR(sigma, i, j) = -MATRIX_AT_PTR(sigma, i, j);
        }
    }
    //  ---- Free alocated memory ----
    freeMatrix(&temp);
    freeMatrix(&diagonalVotesPerBallot);
    free(groupVotesPerBallot);
    // --- ... --- //
}

/**
 * @brief Computes the parameters of the unconditional probability WITH the last candidate
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution.
 *
 * @param[in] b The index of the ballot box
 * @param[in] g The index of the group
 * @param[in] *probabilitiesNotReduced Matrix of dimension (gxc) with the probabilities of each group and candidate,
 * @param[in, out] *mu An array of size c to store the results of the average.
 * @param[in, out] *sigma A matrix of size (c, c) to store the sigma matrix.
 *
 *
 * @return void. Results to be written on mu and sigma.
 *
 */

void getParamsFull(int b, const Matrix *probabilitiesNotReduced, double *mu, Matrix *sigma)
{

    // ---- Check parameters ---- //
    if (probabilitiesNotReduced->cols != TOTAL_CANDIDATES)
    {
        fprintf(stderr, "The probability matrix handed should consider the C candidates, but it has %d columns.\n",
                probabilitiesNotReduced->cols);
        exit(EXIT_FAILURE);
    }
    if (W == NULL && X == NULL)
    {
        fprintf(stderr, "The `w` and `x` matrices aren't defined.\n");
        exit(EXIT_FAILURE);
    }
    // --- ... --- //

    // --- Calculations --- //
    // ---- The votes that a group has made on a given ballot ----
    double *groupVotesPerBallot = getRow(W, b);

    // ---- Computation of mu ----
    // ---- Performing the matrix multiplication of p^T * w_b
    cblas_dgemv(CblasRowMajor,                 // Row-major order
                CblasTrans,                    // Transpose the matrix (p^T)
                TOTAL_GROUPS,                  // Rows of original matrix (G)
                TOTAL_CANDIDATES,              // Columns of original matrix (C)
                1.0,                           // alpha
                probabilitiesNotReduced->data, // Matrix p
                TOTAL_CANDIDATES,              // Leading dimension (rows in original p)
                groupVotesPerBallot,           // Vector w_b
                1,                             // Increment for x
                0.0,                           // beta
                mu,                            // Output vector μ
                1                              // Increment for y
    );

    // ---- Computation of sigma ----
    // ---- Get a diagonal matrix with the group votes on a given ballot ----
    Matrix diagonalVotesPerBallot = createDiagonalMatrix(groupVotesPerBallot, TOTAL_GROUPS);
    // ---- Temporary matrix to store results ----
    Matrix temp = createMatrix(TOTAL_CANDIDATES, TOTAL_GROUPS);
    // ---- Calculates the matrix multiplication of p^T * diag(w_b); result must be (C-1 x G) ----
    cblas_dgemm(CblasRowMajor,                 // Row-major order
                CblasTrans,                    // No transpose (Matrix A => p^T)
                CblasNoTrans,                  // No transpose
                TOTAL_CANDIDATES,              // Rows
                TOTAL_GROUPS,                  // Columns
                TOTAL_GROUPS,                  // Inner dimension
                1.0,                           // alpha
                probabilitiesNotReduced->data, // Matrix p
                TOTAL_CANDIDATES,              // Leading dimension
                diagonalVotesPerBallot.data,   // Matrix diag(w_b)
                TOTAL_GROUPS,                  // Leading dimension
                0.0,                           // beta
                temp.data,                     // Output matrix
                TOTAL_GROUPS);                 // Leading dimension

    // ---- Calculates the matrix multiplication of (p^T * diag(w_b)) * p; result must be (C x C) ----
    cblas_dgemm(CblasRowMajor,                 // Row-major order
                CblasNoTrans,                  // No transpose
                CblasNoTrans,                  // No transpose
                TOTAL_CANDIDATES,              // Rows
                TOTAL_CANDIDATES,              // Columns
                TOTAL_GROUPS,                  // Inner dimension
                1.0,                           // alpha
                temp.data,                     // Matrix (p^T * diag(w_b))
                TOTAL_GROUPS,                  // Leading dimension
                probabilitiesNotReduced->data, // Matrix p
                TOTAL_CANDIDATES,              // Leading dimension
                0.0,                           // beta
                sigma->data,                   // Output Σ_b
                TOTAL_CANDIDATES);             // Leading dimension

    // ---- Substract the diagonal with the average ----
    // ---- Note: This could be optimized with a cBLAS call too ----

    for (int j = 0; j < TOTAL_CANDIDATES; j++)
    { // ---- For each candidate
        for (int i = 0; i < TOTAL_CANDIDATES; i++)
        { // ---- For each candidate given another candidate
            if (i == j)
            { // ---- If it corresponds to a diagonal, substract diagonal
                MATRIX_AT_PTR(sigma, i, j) = mu[i] - MATRIX_AT_PTR(sigma, i, j);
                continue;
            }
            MATRIX_AT_PTR(sigma, i, j) = -MATRIX_AT_PTR(sigma, i, j);
        }
    }
    //  ---- Free alocated memory ----
    freeMatrix(&temp);
    freeMatrix(&diagonalVotesPerBallot);
    free(groupVotesPerBallot);
    // --- ... --- //
}

/**
 * @brief Computes the parameters of the conditional probability WITHOUT the last candidate
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution conditional to the results
 * of a group.
 *
 * @param[in] b The index of the ballot box
 * @param[in] g The index of the group
 * @param[in] *probabilitiesReduce Matrix of dimension (gxc-1) with the probabilities of each group and candidate,
 * except the last one. Consider that the information of the last candidate is redundant. The probability matrix could
 * be reduced with the "removeRows()" function from matrixUtils.
 * @param[in, out] *newMu An array of size c-1 to store the results of the average.
 * @param[in, out] *newSigma A matrix of size (c-1, c-1) to store the sigma matrix.
 *
 * @warning Remember to eliminate one dimension for candidates.
 *
 * @return void. Results to be written on mu and sigma.
 *
 */

void getAverageConditional(int b, const Matrix *probabilitiesReduced, Matrix *conditionalMu, Matrix **conditionalSigma)
{
    // ---- Get the parameters of the unconditional probability ---- //
    double *newMu = (double *)malloc((TOTAL_CANDIDATES - 1) * sizeof(double));
    Matrix newSigma = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1);
    getParams(b, probabilitiesReduced, newMu, &newSigma);
    // ---- ... ----

    // ---- Computation for mu ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
        { // ---- For each candidate given a group
            MATRIX_AT_PTR(conditionalMu, g, c) = newMu[c] - MATRIX_AT_PTR(probabilitiesReduced, g, c);
        }
    }
    // ---- The original mu isn't needed anymore ---- //
    free(newMu);
    //  ---- ... ---- //

    // ---- Get the parameters for the conditional sigma ---- //

    // ---- Get the diagonal probabilities ----
    // ---- Create an array of size `TOTAL_GROUPS` that will store the probabilities for a given group ----
    double **probabilitiesForG = (double **)malloc(TOTAL_GROUPS * sizeof(double *));
    // ---- Create an array of size `TOTAL_GROUPS` that will store diagonal matrices with the probabilities ----
    Matrix *diagonalProbabilities = (Matrix *)malloc((TOTAL_GROUPS) * sizeof(Matrix));
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        probabilitiesForG[g] = getRow(probabilitiesReduced, g);
        diagonalProbabilities[g] = createDiagonalMatrix(probabilitiesForG[g], TOTAL_CANDIDATES - 1);
    }
    // --- ... --- //

    // ---- Get the matrix multiplications ---- //
    // ---- This multiplications are esentially outer products ----
    // ---- Create an array of size `TOTAL_GROUPS` that will store each outer product ----
    Matrix *matrixMultiplications = (Matrix *)malloc((TOTAL_GROUPS) * sizeof(Matrix));
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        // ---- Do the outer product and store it in the array ----
        Matrix mult = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1);
        cblas_dger(CblasRowMajor, TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1, 1.0, probabilitiesForG[g], 1,
                   probabilitiesForG[g], 1, mult.data, TOTAL_CANDIDATES - 1);
        matrixMultiplications[g] = mult;
    }
    // --- ... --- //

    // ---- Add the results to the final array of matrices ----
    // ---- Esentially computes: $$\sigma_b = diag(p_{g}^{t})-p^{t}_{g}p_{g}$$ ----
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        for (uint16_t i = 0; i < TOTAL_CANDIDATES - 1; i++)
        { // ---- For each candidate given a group
            for (uint16_t j = 0; j < TOTAL_CANDIDATES - 1; j++)
            { // ---- For each candidate given a group and a candidate
                // ---- Add the multiplication of probabilities ----
                MATRIX_AT_PTR(conditionalSigma[g], i, j) =
                    MATRIX_AT(newSigma, i, j) + MATRIX_AT(matrixMultiplications[g], i, j);
                if (i == j)
                { // ---- If it's a diagonal
                    // ---- Substract the diagonal probabilities ----
                    MATRIX_AT_PTR(conditionalSigma[g], i, j) -= MATRIX_AT(diagonalProbabilities[g], i, j);
                }
            }
        }
        // ---- Free unnecesary space ----
        freeMatrix(&matrixMultiplications[g]);
        freeMatrix(&diagonalProbabilities[g]);
    }
    // ---- Free space ----
    free(matrixMultiplications);
    free(probabilitiesForG);
    freeMatrix(&newSigma);

    // --- ... --- //
}

/**
 * @brief Computes the parameters of the conditional probability WITH the last candidate
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution conditional to the results
 * of a group.
 *
 * @param[in] b The index of the ballot box
 * @param[in] g The index of the group
 * @param[in] *probabilitiesNotReduced Matrix of dimension (gxc) with the probabilities of each group and candidate,
 * @param[in, out] *newMu An array of size c to store the results of the average.
 * @param[in, out] *newSigma A matrix of size (c, c) to store the sigma matrix.
 *
 * @return void. Results to be written on mu and sigma.
 *
 */

void getAverageConditionalFull(int b, const Matrix *probabilitiesNotReduced, Matrix *conditionalMu,
                               Matrix **conditionalSigma)
{
    // ---- Get the parameters of the unconditional probability ---- //
    double *newMu = (double *)malloc((TOTAL_CANDIDATES) * sizeof(double));
    Matrix newSigma = createMatrix(TOTAL_CANDIDATES, TOTAL_CANDIDATES);
    getParamsFull(b, probabilitiesNotReduced, newMu, &newSigma);
    // ---- ... ----

    // ---- Computation for mu ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        { // ---- For each candidate given a group
            MATRIX_AT_PTR(conditionalMu, g, c) = newMu[c] - MATRIX_AT_PTR(probabilitiesNotReduced, g, c);
        }
    }
    // ---- The original mu isn't needed anymore ---- //
    free(newMu);
    //  ---- ... ---- //

    // ---- Get the parameters for the conditional sigma ---- //

    // ---- Get the diagonal probabilities ----
    // ---- Create an array of size `TOTAL_GROUPS` that will store the probabilities for a given group ----
    double **probabilitiesForG = (double **)malloc(TOTAL_GROUPS * sizeof(double *));
    // ---- Create an array of size `TOTAL_GROUPS` that will store diagonal matrices with the probabilities ----
    Matrix *diagonalProbabilities = (Matrix *)malloc((TOTAL_GROUPS) * sizeof(Matrix));
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        probabilitiesForG[g] = getRow(probabilitiesNotReduced, g);
        diagonalProbabilities[g] = createDiagonalMatrix(probabilitiesForG[g], TOTAL_CANDIDATES);
    }
    // --- ... --- //

    // ---- Get the matrix multiplications ---- //
    // ---- This multiplications are esentially outer products ----
    // ---- Create an array of size `TOTAL_GROUPS` that will store each outer product ----
    Matrix *matrixMultiplications = (Matrix *)malloc((TOTAL_GROUPS) * sizeof(Matrix));
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        // ---- Do the outer product and store it in the array ----
        Matrix mult = createMatrix(TOTAL_CANDIDATES, TOTAL_CANDIDATES);
        cblas_dger(CblasRowMajor, TOTAL_CANDIDATES, TOTAL_CANDIDATES, 1.0, probabilitiesForG[g], 1,
                   probabilitiesForG[g], 1, mult.data, TOTAL_CANDIDATES);
        matrixMultiplications[g] = mult;
    }
    // --- ... --- //

    // ---- Add the results to the final array of matrices ----
    // ---- Esentially computes: $$\sigma_b = diag(p_{g}^{t})-p^{t}_{g}p_{g}$$ ----
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
        { // ---- For each candidate given a group
            for (uint16_t j = 0; j < TOTAL_CANDIDATES; j++)
            { // ---- For each candidate given a group and a candidate
                // ---- Add the multiplication of probabilities ----
                MATRIX_AT_PTR(conditionalSigma[g], i, j) =
                    MATRIX_AT(newSigma, i, j) + MATRIX_AT(matrixMultiplications[g], i, j);
                if (i == j)
                { // ---- If it's a diagonal
                    // ---- Substract the diagonal probabilities ----
                    MATRIX_AT_PTR(conditionalSigma[g], i, j) -= MATRIX_AT(diagonalProbabilities[g], i, j);
                }
            }
        }
        // ---- Free unnecesary space ----
        freeMatrix(&matrixMultiplications[g]);
        freeMatrix(&diagonalProbabilities[g]);
    }
    // ---- Free space ----
    free(matrixMultiplications);
    free(probabilitiesForG);
    freeMatrix(&newSigma);

    // --- ... --- //
}

/**
 * @brief Computes the Mahalanobis distance with last candidate adjustment.
 *
 * @param[in] x Pointer to the input feature vector (size C).
 * @param[in] mu Pointer to the mean vector (size C-1).
 * @param[in] inverseSigma Pointer to the inverse covariance matrix (size (C-1) x (C-1)).
 * @param[out] maha Pointer to the resulting Mahalanobis distances (size C).
 * @param[in] size Size of the truncated candidate space (C-1).
 */
void mahanalobis(double *x, double *mu, Matrix *inverseSigma, double *maha, int size)
{
    // ---- Initialize temporary arrays ----
    double diff[size];
    double temp[size];
    double invs_devs[size];

    // ---- Compute the difference vector ---- //
    for (int i = 0; i < size; i++)
    { // ---- For each truncated element
        diff[i] = x[i] - mu[i];
    }
    // --- ... --- //

    // ---- Compute the multiplication  ---- //
    // ---- inverseSigma * diff

    // ---- Note: The upper triangle is filled on the inverse sigma aswell for the Cholensky method.
    cblas_dsymv(CblasRowMajor, CblasLower, size, 1.0, inverseSigma->data, size, diff, 1, 0.0, temp, 1);
    // --- ... --- //

    // ---- Compute Mahalanobis distance (truncated) ---- //
    double mahanobisTruncated = 0.0;
    for (int i = 0; i < size; i++)
    { // ---- For each truncated element
        // ---- The first parenthesis ----
        mahanobisTruncated += diff[i] * temp[i];
        invs_devs[i] = temp[i]; // Store intermediate results
    }
    // --- ... ---//

    // ---- Compute the Mahanalobis distance with the last candidate ---- //
    /*
     * The updated mahanalobis distance for "C" candidates can be written as:
     *
     * $$D_{i}^{2}=D^{2}_{baseline}-\sigma^{-1}(x_i-\mu_i)+diag(\sigma^{-1})$$
     *
     * The baseline would be the mahanobis distance for the "C-1" candidate (mahanobisTruncated)
     * */

    maha[size] = mahanobisTruncated; // Last element is used as a reference
    for (int c = 0; c < size; c++)
    { // ---- For each candidate (considers the last one too)
        maha[c] = mahanobisTruncated - 2 * temp[c] + MATRIX_AT_PTR(inverseSigma, c, c);
    }
    // ---...--- //
}

Matrix getBallotPDF(int b, const Matrix *probabilitiesReduced)
{
    // --- Get the mu and sigma --- //
    Matrix muR = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES - 1);
    Matrix **sigma = (Matrix **)malloc(TOTAL_GROUPS * sizeof(Matrix *));

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        sigma[g] = (Matrix *)malloc(sizeof(Matrix));
        *sigma[g] = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1); // Initialize
    }

    getAverageConditional(b, probabilitiesReduced, &muR, sigma);
    // ---- ... ----

    // ---- Get the inverse matrix for each sigma ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        // ---- Calculates the Cholensky matrix ---- //
        inverseSymmetricPositiveMatrix(sigma[g]);
    }
    // ---- ... ----

    // --- Calculate the mahanalobis distance --- //
    double **mahanalobisDistances = (double **)malloc(TOTAL_GROUPS * sizeof(double *));

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        mahanalobisDistances[g] = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));
        // ---- Get the feature vector (the candidate results) ----
        double *feature = getColumn(X, b);
        // ---- Get the average values for the candidate results ----
        double *muG = getRow(&muR, g);
        // ---- Call the mahanalobis function ----
        mahanalobis(feature, muG, sigma[g], mahanalobisDistances[g], TOTAL_CANDIDATES - 1);
        // ---- Free allocated and temporary values ----
        freeMatrix(sigma[g]);
        free(feature);
        free(muG);
    }
    free(sigma);
    freeMatrix(&muR);
    // --- .... --- //

    // --- Calculate the returning values --- //
    // ---- Create the matrix to return ----
    Matrix toReturn = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        // ---- Initialize variables ----
        double *QC = (double *)calloc(TOTAL_CANDIDATES, sizeof(double)); // Value of Q on candidate C
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        { // ---- For each candidate given a group
            MATRIX_AT(toReturn, g, c) = exp(-0.5 * mahanalobisDistances[g][c]);
        }

        // ---- Free allocated memory ----
        free(mahanalobisDistances[g]);
        free(QC);
    }
    free(mahanalobisDistances); // Might have to remove

    return toReturn;
    // --- ... --- //
}
