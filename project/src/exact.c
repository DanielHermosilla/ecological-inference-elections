#include "matrixUtils.h"
#include "memoizationUtil2.h"
#include <cblas.h>
#include <cstdlib>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf_comb.h>
#include <math.h>
#include <stdint.h>

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

/**
 * @brief Calculate the chained product between probabilities as defined in `a`:
 *
 * Given an `H` element and the `f` index it computes the chained product.
 *
 * @param[in] *hElement A pointer towards the element of `H`. This element is an array with the possible combinations
 * @param[in] *probabilities A pointer toward the probabilities Matrix.
 * @param[in] f The index of `f`
 *
 * @return double: The result of the product
 *
 */
double prod(const size_t *hElement, const Matrix *probabilities, const int f)
{
    double result = 1;
    for (int c = 0; c < TOTAL_CANDIDATES; c++)
    {
        result *= pow(MATRIX_AT_PTR(probabilities, f, c), hElement[c]);
    }
    return result;
}

/**
 * @brief Calculate `a` given the multinomial coefficient, subject to the `H` set:
 *
 * Given an `H` element, it computes its multinomial coefficient subject to the total amount of votes a group received.
 * It yields the following calculation:
 * $$\binom{w_{bf}}{h_1,\dots,h_C}=\frac{w_{bf}!}{h_1!\cdots h_{C}!}$$
 *
 * @param [in] b The index of the ballot.
 * @param [in] f The index of the group
 * @param[in] *hElement A pointer towards the element of `H`. This element is an array with the possible combinations
 *
 * @return double: The result of the multinomial coefficient
 *
 */
double multinomialCoeff(const int b, const int f, const size_t *hElement)
{
    double result = gsl_sf_fact(MATRIX_AT_PTR(W, b, f)); // Calculate w_bf!
    for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
    {
        result /= gsl_sf_fact(hElement[i]); // Divide by each h_i!
    }

    return result;
}

bool ifAllElements(const size_t *hElement, const size_t *kElement)
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

double computeA(const int b, const int f, const size_t *hElement, const Matrix *probabilities)
{
    return multinomialCoeff(b, f, hElement) * prod(hElement, probabilities, f);
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

bool filterCombinations(const size_t *hVector, const int b)
{
    for (int c = 0; c < TOTAL_CANDIDATES; c++)
    {
        if (hVector[c] > (size_t)MATRIX_AT_PTR(X, b, c))
        {
            return false;
        }
    }
    return true;
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

            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                restriction[*size][c] = elementArray[c]; // Copy the element into the array
            }
        }

    } while (gsl_combination_next(noRestriction) == GSL_SUCCESS);

    gsl_combination_free(noRestriction);

    return restriction;
}

void vectorDiff(const size_t *K, const size_t *H, size_t *arr)
{
    for (int c = 0; c < TOTAL_CANDIDATES; c++)
    {
        arr[c] = K[c] - H[c];
    }
}

double recursion2(const int b, const int f, const int g, const int c, size_t *vector, MemoizationTable *memo,
                  const Matrix *probabilities)
{
    // If it's already defined, return the defined value instead; at each iteration this condition will be met more
    // times.
    double currentVal = getMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES);
    if (currentVal != INVALID)
    {
        return currentVal;
    }

    // ---- Base case ----
    if (f == 0)
    {
        for (int k = 0; k < TOTAL_CANDIDATES; k++)
        {
            if (vector[k] != 0)
            {
                setMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES, 0.0);
                return 0.0;
            }
        }
        setMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES, 1.0);
        return 1.0;
    }
    // ---- End base case ----

    // ---- Recursion ----
    // ---- Recursion: definition of H ----
    gsl_combination *H = getH(b, f);
    // ---- Recursion: definition of the resulting variable ----
    double result = 0;
    // ---- Recursion: loop over H, it must be a summatory ----
    do
    {
        // ---- Recursion: get a pointer of an element from H ----
        const size_t *hElement = gsl_combination_data(H);
        // ---- Recursion: handle case when substraction could be negative (*) ----
        if (!ifAllElements(hElement, vector))
            continue;
        // ---- Recursion: compute the value of `a` ----
        double a = computeA(b, f, hElement, probabilities);
        // ---- Recursion: define the vector k-h ----
        size_t *substractionVector = malloc(TOTAL_CANDIDATES * sizeof(int));
        vectorDiff(vector, hElement, substractionVector);

        // ---- Recursion: handle case of the additional multiplication ----
        if (f == g)
        {
            double division = MATRIX_AT_PTR(probabilities, f, c) * MATRIX_AT_PTR(W, b, f);
            result += (recursion2(b, f - 1, g, c, substractionVector, memo, probabilities) * a) / division;
        }
        else
        {
            result += recursion2(b, f - 1, g, c, substractionVector, memo, probabilities) * a;
        }
        // ---- Recursion: free the substracting vector, it won't be used again ----
        free(substractionVector);
        // ---- Note: the element is defined as a summatory from all elements from H_bf that accomplish (*) ----

    } while (gsl_combination_next(H) == GSL_SUCCESS); // Recursion: condition to loop over each H element

    // ---- Recursion: final step, set the hash table to the final value and return it
    setMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES, result);
    return result;
    // ---- End recursion ----
}

void getUrecursive2(MemoizationTable *memo, const Matrix *probabilities, int g, int c)
{
    /* Here is where the loop starts and the values are defined */
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (int f = 0; f < TOTAL_GROUPS; f++)
        {
            size_t sizeK;
            size_t **K = getK(b, f, &sizeK); // K_{bf}

            double division = MATRIX_AT_PTR(probabilities, f, c) *
                              MATRIX_AT_PTR(W, b, f); // Computes now the division to avoid multiple divisions.

            // ---- For each element of K_{bf} ----
            for (int k = 0; k < sizeK; k++) // For each element of K_{bf}
            {                               // Loop over K
                size_t *kElement = K[k];
                gsl_combination *H = getH(b, f);

                // ---- Base case ----
                if (f == 0)
                {
                    // If `k` is a null vector, then return "1", else, "0"
                    bool passed = false;
                    for (int k = 0; k < TOTAL_CANDIDATES; k++)
                    {
                        if (kElement[k] != 0 && !passed)
                        {
                            setMemoValue(memo, b, f, g, c, kElement, TOTAL_CANDIDATES, 0.0);
                            passed = true;
                        }
                    }
                    if (!passed)
                        setMemoValue(memo, b, f, g, c, kElement, TOTAL_CANDIDATES, 1.0);
                }
                // ---- End base case ----
                else
                {

                    // ---- For each element of H_{bf} ----
                    // The gsl_combination library offers a way to iterate over the elements of a set, hence, the
                    // iteration with H will be made with that method in mind on a do-while loop.
                    do
                    {
                        const size_t *hElement = gsl_combination_data(H);
                        if (!ifAllElements(hElement,
                                           kElement)) // Handle case when substraction vector could be negative

                            continue;
                        double a = computeA(b, f, hElement, probabilities);
                        size_t *substractionVector = malloc(TOTAL_CANDIDATES * sizeof(int));
                        vectorDiff(kElement, hElement, substractionVector);

                        double value = getMemoValue(memo, b, f, g, c, kElement, TOTAL_CANDIDATES);
                        if (value == INVALID) // If the entry doesn't exist yet
                        {
                            value = 0.0;
                        }
                        if (f == g)
                        {
                            double division = MATRIX_AT_PTR(probabilities, f, c) * MATRIX_AT_PTR(W, b, f);
                            value += getMemoValue(memo, b, f - 1, g, c, substractionVector, TOTAL_CANDIDATES) * a *
                                     hElement[c] / division;
                        }
                        else
                        {
                            value += getMemoValue(memo, b, f - 1, g, c, substractionVector, TOTAL_CANDIDATES) * a;
                        }
                        setMemoValue(memo, b, f, g, c, kElement, TOTAL_CANDIDATES, value);
                        free(substractionVector);

                    } while (gsl_combination_next(H) == GSL_SUCCESS); // Loop over each H element
                }
            }
        }
    }
}
void computeQExact()
{

    MemoizationTable *table = initMemo();
}
