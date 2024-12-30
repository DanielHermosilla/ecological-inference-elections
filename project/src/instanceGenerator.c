#include "matrixUtils.h"
#include <stdint.h>
#include <stdlib.h>

// Helper function to get a valid integer from the user
int getInput(const char *prompt)
{
    int value;
    int valid;

    while (1)
    { // Loop until valid input
        printf("%s", prompt);
        valid = scanf("%d", &value);
        if (valid == 1)
        {
            break; // Valid integer entered, exit the loop
        }
        else
        {
            printf("Invalid input. Answer must be an integer.\n");
            while (getchar() != '\n')
                ; // This would clear word by word the input buffer until it's empty.
        }
    }
    return value;
}

/**
 * @brief Generates an instance according some parameters
 *
 * Given some input parameters, it generates random matrices that represent a votation
 *
 * @return TODO; maybe an array that points to the matrices.
 *
 * @note This would be for debugging,
 *
 * @see getInitialP() for getting initial probabilities. Group proportional method is recommended.
 *
 * @warning
 * - Pointers shouldn't be NULL.
 * - `x` and `w` dimensions must be coherent.
 *
 */

void createInstance()
{

    printf("Welcome to the instance maker!\n");

    int totalvotes = getInput("Please enter the total amount of votes: \n");
    int totalcandidates = getInput("Please enter the total amount of candidates: \n");
    int totalgroups = getInput("Please enter the total amount of demographic groups: \n");
    int totalballots = getInput("Please enter the total amount of ballot boxes: \n");

    printf("\nInstance Summary:\n");
    printf("Total Votes: %d\n", totalvotes);
    printf("Total Candidates: %d\n", totalcandidates);
    printf("Total Groups: %d\n", totalgroups);
    printf("Total Ballots: %d\n", totalballots);
}

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

void distributeTotalRandomly(int total, int size, double *array)
{
    for (int i = 0; i < total; i++)
    { // The loop will assign a vote to a random index

        int index = rand() % size; // Choose a random index
        array[index]++;            // Adds a vote to that random index
    }
}

void generateVotes(const int *totalvotes, const int *totalcandidates, const int *totalgroups, const int *totalballots)
{

    Matrix x = createMatrix(*totalcandidates, *totalballots); // The dimensions are (cxb)
    Matrix w = createMatrix(*totalballots, *totalgroups);     // The dimensions are (bxg)

    // The sum of columns in x MUST be the same as the sum of rows in w.

    // Firstly, the total amount of votes per each ballot should be random.
    double votesPerBallot[*totalballots];
    makeArray(votesPerBallot, *totalballots, 0.0);
    distributeTotalRandomly(*totalvotes, *totalballots, votesPerBallot);

    // Now, votesPerBallot generated a random total amount of votes per ballot.
    // Filling matrix x:
#pragma omp parallel for
    for (int b = 0; b < *totalballots; b++)
    {
        double candidateVotes[*totalcandidates];
        double groupVotes[*totalgroups];
        distributeTotalRandomly(
            votesPerBallot[b], *totalcandidates,
            candidateVotes); // Distributes the total amount of votes of a candidate randomly per ballot box.
        distributeTotalRandomly(
            votesPerBallot[b], *totalgroups,
            groupVotes); // Distributes the total amount of votes of a group randomly perr ballot box.
        for (int c = 0; c < *totalcandidates; c++)
        {
            MATRIX_AT(x, c, b) = candidateVotes[c];
        }
        for (int g = 0; g < *totalgroups; g++)
        {
            MATRIX_AT(w, b, g) = groupVotes[g];
        }
    }
}

int main()
{
    return 1;
}

// Generate a random number in the range [min, max]
// int rd_num = rand_r(&seed) % (max - min + 1) + min;
// printf("%d ", rd_num);
