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
bool filterCombinations(const size_t *element, const uint32_t b)
{

    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    {
        if (element[c] < MATRIX_AT_PTR(X, c, b))
        {
            return false;
        }
    }
    return true;
}
void defineK()
{
    Matrix K = createMatrix(W->rows, W->cols);
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            if (g == 0)
            {
                MATRIX_AT(K, b, g) = MATRIX_AT_PTR(W, b, g);
            }
            MATRIX_AT(K, b, g) += MATRIX_AT(K, b, g - 1);
        }
    }
}

gsl_combination *getH(int b, int f)
{
    gsl_combination *result = gsl_combination_calloc(MATRIX_AT_PTR(W, b, f), TOTAL_CANDIDATES);
    return &result;
}

gsl_combination *getK(int b, int f, size_t *size)
{

    if (f > TOTAL_GROUPS)
    {
        printf("An incorrect amount of groups was passed to the `K` set");
        exit(EXIT_FAILURE);
    }

    size_t **restriction = NULL;
    *size = 0;

    int totalVotes = 0;
    for (int i = 0; i < f; i++)
    {
        totalVotes += MATRIX_AT_PTR(W, b, i);
    }
    gsl_combination *noRestriction = getH(totalVotes, TOTAL_CANDIDATES);

    do
    {
        const size_t *elementArray = gsl_combination_data(noRestriction);

        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {
        }

    } while (gsl_combination_next(noRestriction) == GSL_SUCCESS);
    return &result;
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
