#include "matrixUtils.h"
#include <cblas.h>
#include <cstdlib>
#include <gsl/gsl_combination.h>

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

Matrix *KSET = NULL;

typedef struct
{
    size_t n;
    size_t k;
    size_t *data;
} gsl_combination;
// Funciones utiles;
// gsl_sl_fact(unsigned int n) Factorial
//
// Combinatoria: Está en el siguiente structure
// typedef struct
// {
// size_t n;
// size_t k;
// size_t *data;
//} gsl_combination;

/*
gsl_combination *gsl_combination_calloc(size_t n, size_t k)

    This function allocates memory for a new combination with parameters n, k and initializes it to the
lexicographically first combination. A null pointer is returned if insufficient memory is available to create the
combination.
*/

/*
void gsl_combination_free(gsl_combination *c)

    This function frees all the memory used by the combination c.
*/
bool ifAllElements(size_t *hElement, size_t *kElement)
{
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    {
        if (hElement[c] > kElement[c])
        {
            return false;
        }
    }
    return true;
}

/**
 * @brief Calculate a given `H` set.
 *
 * Given the index of the set, it returns the `H` set as a gsl_combination struct. Note that it could be used with the
 * gsl_combination functions since it doesn't have restrictions.
 *
 * @param[in] b Index that represents the corresponding ballot.
 * @param[in] f Index that represents the f group (starting at 0).
 *
 * @return gsl_combination *: A pointer towards the struct of a gsl_combination.
 *
 * @warning
 * - The `f` index cannot be bigger than the total amount of groups.
 * - The `b` index cannot be bigger than the total ballots.
 *
 */
gsl_combination *getH(int b, int f)
{
    gsl_combination *result = gsl_combination_calloc(MATRIX_AT_PTR(W, b, f), TOTAL_CANDIDATES);
    return result;
}

/**
 * @brief Calculate a given `K` set.
 *
 * Given the index of the set, it returns the `K` set and its size. Note that it will be an array of arrays. Both arrays
 * are size_t.
 *
 * @param[in] b Index that represents the corresponding ballot.
 * @param[in] f Index that represents the first `f` groups (starting at 0).
 * @param[in, out] size The size of the new array.
 *
 * @return size_t **: A pointer towards the elements of the array of arrays.
 *
 * @note TODO: The while loop could be optimized.
 *
 * @warning
 * - The `f` index cannot be bigger than the total amount of groups.
 * - The `b` index cannot be bigger than the total ballots.
 *
 */

size_t **getK(const int b, const int f, size_t *size) // combination[i][j]; i the combination, j the element
{
    if (f >= TOTAL_GROUPS || b >= TOTAL_BALLOTS)
    {
        printf("An incorrect index was passed to the `K` set");
        exit(EXIT_FAILURE);
    }

    size_t **restriction = NULL;
    *size = 0;

    int totalVotes = 0; // The total votes accumulated for the first `f` groups
    for (int i = 0; i < f; i++)
    {
        totalVotes += MATRIX_AT_PTR(W, b, i);
    }
    gsl_combination *noRestriction = getH(totalVotes, TOTAL_CANDIDATES);

    do // It could be optimized by stepping on a big amount of elements.
    {
        const size_t *elementArray = gsl_combination_data(noRestriction);

        if (filterCombinations(elementArray, b))
        {
            restriction = realloc(restriction, (*size + 1) * sizeof(size_t *)); // The size of the whole array
            restriction[*size] = malloc(TOTAL_CANDIDATES * sizeof(size_t)); // The size of a single element (constant)
        }

        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            restriction[*size][c] = elementArray[c]; // Copy the element into the array
        }

    } while (gsl_combination_next(noRestriction) == GSL_SUCCESS);

    gsl_combination_free(noRestriction);

    return restriction;
}

double prod(const int i, const int n, const double *array)
{
    double result = 1;
    for (int k = i; k < n; k++)
    {
        result *= array[k];
    }
    return result;
}

bool checkCondition(const int b, const int f)
{
    gsl_combination *currentH = getH(b, f);
}

void computeExact()
{

    int sizeK = W->rows * W->cols;
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (int f = 0; f < TOTAL_GROUPS; f++)
        {
            for (int k = 0; k < sizeK; k++)
            { // That's the size of K, this should be a set loop
                for (int h = 0; h < sizeK; h++)
                {
                    if (MATRIX_AT(H, b, g))
                }
            }
        }
    }
}
