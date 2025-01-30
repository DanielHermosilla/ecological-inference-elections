#include "matrixUtils.h"
#include <ctype.h> // For tolower()
#include <stdint.h>
#include <stdlib.h>

#ifndef INSTANCE_H
#define INSTANCE_H

/**
 * @brief Helper function to get a valid integer from the user
 *
 * Ensures that the user only gives an integer.
 *
 * @param[in] *prompt The prompt for the input
 *
 * @return the value that was written by the user.
 */

int getInput(const char *prompt);

/**
 * @brief Distributes, randomly, amount of votes into an array.
 *
 * Given a total amount of votes, it fills an array with those votes randomly.
 *
 * @param[in] total The total amount of votes
 * @param[in] size The size of the array
 * @param[in, out] *array The array to write the results
 *
 * @return void
 */

void distributeTotalRandomly(int total, int size, double *array);

/**
 * @brief Helper function for generating the random instances of voting.
 *
 * The function will be in charge of filling each ballot randomly while making sure that the sums of `x` columns
 * must be the same as the sum of `w` rows.
 *
 * @param[in, out] x A pointer to the matrix `x`
 * @param[in, out] w A pointer to the matrix `w`
 * @param[in] totalvotes The total amount of votes
 * @param[in] totalcandidates The total amount of candidates
 * @param[in] totalgroups The total amount of groups
 * @param[in] totalballots The total amount of ballots
 *
 * @return void
 */

void generateVotes(Matrix *x, Matrix *w, const int *totalvotes, const int *totalcandidates, const int *totalgroups,
                   const int *totalballots);

/**
 * @brief Generates an instance according some parameters
 *
 * Given some input parameters, it generates random matrices that represent a votation
 *
 * @param[in, out] x A pointer to the matrix `x` to be generated. I would sugguest to just pass an empty matrix pointer.
 * @param[in, out] w A pointer to the matrix `w` to be generated. I would sugguest to just pass an empty matrix pointer.
 * @param[in] seed The seed for making reproducible results.
 * @param[in] method The method for generating votes, it could either be "uniform" or "multinomial".
 *
 * @return void. Results written at the matrix pointers
 *
 * @note This would be mainly for debugging,
 */

void createInstance(Matrix *x, Matrix *w, const int seed, const char method);

#endif // INSTANCE_H
