// El header del paquete, acá llamamos las funciones que pasarán a R
#ifndef MAIN_H
#define MAIN_H

#ifdef __cplusplus

extern "C"
{
#endif

#include "exact.h"
#include "globals.h"
#include "hitAndRun.h"
// #include "instanceGenerator.h"
#include "multinomial.h"
#include "utils_matrix.h"
// #include "multivariate-cdf.h"
#include "multivariate-pdf.h"
#include "utils_file.h"

    // ---- Define the structure to store the function pointer ---- //
    // ----  It is defined here since it's not used globally ----
    typedef struct
    {
        double *(*computeQ)(const Matrix *, QMethodInput); // Function pointer for computing q
        QMethodInput params;                               // Holds method-specific parameters
    } QMethodConfig;
    // ---...--- //

    /**
     * @brief Yields the global parameters of the process. Usually this should be done once for avoiding
     * computing a loop over ballots. It also changes the parameters in case it's called with other `x` and `w` matrix.
     *
     * Gets the total amount of votes in the process. This should only be donce once for avoiding computing loops
     * over ballots.
     *
     * @param[in] x Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
     * @param[in] w Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
     *
     * @return void. Will edit the static values in the file
     *
     * @note This should only be used once, later to be declared as a static value in the program.
     *
     * @warning
     * - Pointers shouldn't be NULL.
     * - `x` and `w` dimensions must be coherent.
     *
     */
    void setParameters(Matrix *x, Matrix *w);

    /**
     * @brief Computes the initial probability of the EM algoritm.
     *
     * Given the observables results, it computes a convenient initial "p" value for initiating the
     * algorithm. Currently it supports the "uniform", "group_proportional" and "proportional" methods.
     *
     * @param[in] p_method The method for calculating the initial parameter. Currently it supports "uniform",
     * "group_proportional" and "proportional" methods.
     *
     * @return Matrix of dimension (gxc) with the initial probability for each demographic group "g" voting for a given
     * candidate "c".
     * @note This should be used only that the first iteration of the EM-algorithm.
     * @warning
     * - Pointers shouldn't be NULL.
     * - `x` and `w` dimensions must be coherent.
     *
     */
    Matrix getInitialP(const char *p_method);

    /*
     * @brief Computes the optimal solution for the `M` step
     *
     * Given the conditional probability and the votations per demographic group, it calculates the new probability for
     * the next iteration.
     *
     * @param[in] q Array of matrices of dimension (bxgxc) that represents the probability that a voter of group "g" in
     * ballot box "b" voted for candidate "c" conditional on the observed result.
     *
     * @return A matrix with the optimal probabilities according maximizing the Log-likelihood.
     *
     * @see getInitialP() for getting initial probabilities. This method is recommended to be used exclusively for the
     * EM Algorithm, unless there's a starting "q" to start with.
     *
     */

    //  Matrix getP(const double *q);

    /**
     * @brief Implements the whole EM algorithm.
     *
     * Given a method for estimating "q", it calculates the EM until it converges to arbitrary parameters. As of in the
     * paper, it currently supports Hit and Run, Multinomial, MVN CDF and MVN PDF methods.
     *
     * @param[in] currentP Matrix of dimension (cxg) with the initial probabilities for the first iteration.
     * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
     * and Run", "Multinomial", "MVN CDF", "MVN PDF" and "Exact" methods.
     * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
     * @param[in] maxIter Integer with a threshold of maximum iterations. Usually it's set to 100.
     * @param[in] verbose Wether to verbose useful outputs.
     *
     * @return Matrix: A matrix with the final probabilities. In case it doesn't converges, it returns the last
     * probability that was computed
     *
     * @note This is the main function that calls every other function for "q"
     *
     * @see getInitialP() for getting initial probabilities. Group proportional method is recommended.
     *
     * @warning
     * - Pointers shouldn't be NULL.
     * - `x` and `w` dimensions must be coherent.
     *
     */
    Matrix EMAlgoritm(Matrix *currentP, const char *q_method, const double convergence, const int maxIter,
                      const bool verbose, double *time, int *iterTotal, double **logLLarr, QMethodInput params);

    /**
     * @brief Checks if a candidate didn't receive any votes.
     *
     * Given an array of size TOTAL_CANDIDATES, it sets to "1" the index where a possible candidate haven't received any
     * vote. It also returns a boolean indicating whether a candidate hasn't receive any vote
     *
     * @param[in,out] *canArray Array of size TOTAL_CANDIDATES full of zeroes, indicating with a "1" on the index where
     * a given candidate haven't received a vote
     *
     * @return bool: A boolean that shows if it exists a candidate with no votes
     *
     */
    //   bool noVotes(int *canArray);
#ifdef __cplusplus
}
#endif
#endif // UTIL_H
