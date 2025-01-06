#include "exact.h"
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

#define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[((bIdx) * (G) * (C)) + ((gIdx) * (C)) + (cIdx)])
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
 * @brief Converts a type `double` array to a `size_t`
 *
 * Given an array that is made with `double` elements, it fills another one with `size_t` elements. Note that the
 * decimal values will be truncated.
 *
 * @param[in] *doubleVector A pointer towards the vector that has the values in `double`.
 * @param[in] *sizeTVector A pointer toward the vector with values of type `size_t`.
 * @param[in] size The size of the array
 *
 * @return void
 *
 */
void convertDoubleToSizeT(const double *doubleVector, size_t *sizeTVector, int size)
{
    for (int i = 0; i < size; i++)
    {
        sizeTVector[i] = (size_t)doubleVector[i];
    }
}

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

/**
 * @brief Checks if all elements of `K` are bigger than the elements from `H` to ensure non negativity.
 *
 * Given vectors `H` and `K` it iterates between all its values to see the sign of the difference
 *
 * @param[in] *hElement A pointer to the array `H`
 * @param[in] *kElement A pointer to the array `K`
 *
 * @return Boolean value that tells if the condition is met.
 *
 */

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

/**
 * @brief Computes the `A` term from the pseudocode
 *
 * Given the vector `H` and the index from the ballot box and cummulative group, it computes the value from `A`, defined
 * as a product of inner products and binomial coefficients.
 *
 * @param[in] b Index that represents the corresponding ballot.
 * @param[in] f Index that represents the f group (starting at 0).
 * @param[in] *hElement Pointer to the `H` vector
 * @param[in] *probabilities Pointer to the probabilities matrix.
 *
 * @return The value of `A`
 *
 */

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

/**
 * @brief Tells if a combination is valid according the `K` definition
 *
 * Given the ballot box identificator, it checks if the condition for joining `K` is met.
 *
 * @param[in] *hVector
 * @param[in] b Index that represents the ballot box.
 *
 * @return A boolean that tells if the condition is met
 *
 *
 */

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

/**
 * @brief Calculate the difference between two vectors.
 *
 * Utility function to get the difference between two 1-dimensional arrays.
 *
 * @param[in] *K A pointer to the first array
 * @param[in] *H A pointer to the second array
 * @param[in, out] *arr The result of the operation
 *
 * @return void
 *
 */

void vectorDiff(const size_t *K, const size_t *H, size_t *arr)
{
    for (int c = 0; c < TOTAL_CANDIDATES; c++)
    { // ---- For each candidate
        arr[c] = K[c] - H[c];
    }
}

/**
 * @brief Calculate the value of `u_{bfgc}(k)`.
 *
 * Given the index from `u` it returns its value. If it hasn't been calculated yet, it starts iterating the recursion
 * for getting a final answer.
 *
 * @param[in] b Index that represents the corresponding ballot.
 * @param[in] f Index that represents the first `f` groups (starting at 0).
 * @param[in] g Index that represents the `g` group (starting at 0).
 * @param[in] c Index that represents the `c` candidate (starting at 0).
 * @param[in] *vector A pointer to the vector that functions as a key.
 * @param[in, out] *memo A pointer to the hash table.
 * @param[in] *probabilities A pointer to the matrix with the probabilities.
 *
 * @return double: The value from `u_{bfgc}(k)`
 *
 * @warning
 * - The `f` index cannot be bigger than the total amount of groups.
 * - The `b` index cannot be bigger than the total ballots.
 *
 */

double recursion(const int b, const int f, const int g, const int c, size_t *vector, MemoizationTable *memo,
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
        { // ---- For each candidate
            if (vector[k] != 0)
            { // ---- If the vector contains an non zero element; return 0.0
                setMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES, 0.0);
                return 0.0;
            }
        }
        // ---- Up to this point, the vector should be full of zeros.
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
            result += (recursion(b, f - 1, g, c, substractionVector, memo, probabilities) * a) / division;
        }
        else
        {
            result += recursion(b, f - 1, g, c, substractionVector, memo, probabilities) * a;
        }
        // ---- Recursion: free the substracting vector, it won't be used again ----
        free(substractionVector);
        // ---- Note: the element is defined as a summatory from all elements from H_bf that accomplish (*) ----

    } while (gsl_combination_next(H) == GSL_SUCCESS); // Recursion: condition to loop over each H element

    // ---- Recursion: free the H set that was created ----
    // TODO: Maybe precompute the set to avoid multiple computations
    gsl_combination_free(H);
    // ---- Recursion: final step, set the hash table to the final value and return it
    setMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES, result);
    // ---- Note: the vector that was used as a key will be free'd with the table, since it's a pointer.
    return result;
    // ---- End recursion ----
}

/**
 * @brief Calculate the value of `q_{bgc}`.
 *
 * It calculates all of the values for q_{bgc} by the definition on the paper. It returns the array of type `double`.
 *
 * @param[in] *probabilities A pointer to the matrix with the probabilities.
 *
 * @return *double: A pointer toward the array.
 *
 * @note: A single pointer is used to store the array continously. This is for using cBLAS operations later.
 *
 */

double *computeQExact(const Matrix *probabilities)
{
    // ---- Initialize the hash table
    MemoizationTable *table = initMemo();
    // ---- Initialize the array to return
    double *array2 = (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double));
    // ---- Initialize the loop for computing each element of q_{bgc}
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For each ballot box
        for (int g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group
            // ---- Initialize an array of arrays to store the candidates votes on `size_t`
            size_t **candidateArrays = malloc(TOTAL_CANDIDATES * sizeof(size_t *));
            // ---- Initialize the denominator variable to avoid multiple computations
            double den = 0;
            for (int c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate
                // ---- Initialize an array for the outer array
                candidateArrays[c] = malloc(TOTAL_CANDIDATES * sizeof(size_t));
                // ---- Convert the array to `size_t`
                convertDoubleToSizeT(&MATRIX_AT_PTR(X, c, b), candidateArrays[c], TOTAL_CANDIDATES);
                // ---- Add the values of the denominator
                den += recursion(b, TOTAL_GROUPS, g, c, candidateArrays[c], table, probabilities) *
                       MATRIX_AT_PTR(probabilities, g, c);
            }
            for (int c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate
                // ---- Compute the numerator; the recursion shouldn't be ran since the denominator initialized the
                // hashed values.
                double num = recursion(b, TOTAL_GROUPS, g, c, candidateArrays[c], table, probabilities) *
                             MATRIX_AT_PTR(probabilities, g, c);
                // ---- Store the resulting values as q_{bgc}
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = num / den;
            }
            // ---- Frees the array of arrays outer pointer. This should be safe since the keys of the Hash Table has a
            // pointer towards the inner arrays. Something that could be done is to store this value since it's going to
            // be reused for each EM iteration.
            free(candidateArrays);
        }
    }
    freeMemo(table);
    return array2;
}
