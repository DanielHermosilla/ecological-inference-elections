#include "matrixUtils.h"
#include <ctype.h> // For tolower()
#include <stdint.h>
#include <stdlib.h>

/**
 * @brief Helper function to get a valid integer from the user
 *
 * Ensures that the user only gives an integer.
 *
 * @param[in] *prompt The prompt for the input
 *
 * @return the value that was written by the user.
 */

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
                   const int *totalballots)
{
    // The sum of columns in x MUST be the same as the sum of rows in w.

    // Firstly, the total amount of votes per each ballot should be random.
    double votesPerBallot[*totalballots];
    makeArray(votesPerBallot, *totalballots, 0.0);
    distributeTotalRandomly(*totalvotes, *totalballots, votesPerBallot);

    // Now, votesPerBallot generated a random total amount of votes per ballot.
    // Filling matrix x:
    // #pragma omp parallel for
    for (int b = 0; b < *totalballots; b++)
    {
        double candidateVotes[*totalcandidates];
        double groupVotes[*totalgroups];
        makeArray(candidateVotes, *totalcandidates, 0.0);
        makeArray(groupVotes, *totalgroups, 0.0);
        distributeTotalRandomly(
            votesPerBallot[b], *totalcandidates,
            candidateVotes); // Distributes the total amount of votes of a candidate randomly per ballot box.
        distributeTotalRandomly(votesPerBallot[b], *totalgroups,
                                groupVotes); // Distributes the total amount of votes of a group randomly perr ballot
                                             // box. #pragma omp critical
        //        {
        // printf("\nFor the ballot %d, the candidate distribution will be:\n", b);
        for (int c = 0; c < *totalcandidates; c++)
        {
            MATRIX_AT_PTR(x, c, b) = candidateVotes[c];
            // printf("%.1f\t", candidateVotes[c]);
        }
        // printf("\nAnd the group distribution will be:\n");
        for (int g = 0; g < *totalgroups; g++)
        {
            MATRIX_AT_PTR(w, b, g) = groupVotes[g];
            // printf("%.1f\t", groupVotes[g]);
        }
        // printf("\nIn theory, both should sum to %.1f\n\n", votesPerBallot[b]);
    }
    //    }
}

/**
 * @brief Generates an instance according some parameters
 *
 * Given some input parameters, it generates random matrices that represent a votation
 *
 * @param[in, out] x A pointer to the matrix `x` to be generated. I would sugguest to just pass an empty matrix pointer.
 * @param[in, out] w A pointer to the matrix `w` to be generated. I would sugguest to just pass an empty matrix pointer.
 * @param[in] seed The seed for making reproducible results.
 *
 * @return void. Results written at the matrix pointers
 *
 * @note This would be mainly for debugging,
 */

void createInstance(Matrix *x, Matrix *w, const int seed)
{
    srand(seed);

    printf("Welcome to the instance maker!\n");

    int totalvotes = getInput("Please enter the total amount of votes: \n");
    int totalcandidates = getInput("Please enter the total amount of candidates: \n");
    int totalgroups = getInput("Please enter the total amount of demographic groups: \n");
    int totalballots = getInput("Please enter the total amount of ballot boxes: \n");

    if (totalgroups >= totalvotes || totalcandidates >= totalvotes)
    {
        fprintf(stderr,
                "There are more groups (or candidates) than the total amount of votes, please eliminate redundancy");
        exit(EXIT_FAILURE);
    }

    // Validate initial state
    if (x->data != NULL || w->data != NULL)
    {
        fprintf(stderr, "Matrices were not empty before initialization.\n");
        exit(EXIT_FAILURE);
    }

    x->data = (double *)malloc(totalcandidates * totalballots * sizeof(double));
    if (!x->data)
    {
        fprintf(stderr, "Failed to allocate matrix X.\n");
        exit(EXIT_FAILURE);
    }
    x->rows = totalcandidates;
    x->cols = totalballots;

    w->data = (double *)malloc(totalgroups * totalballots * sizeof(double));
    if (!w->data)
    {
        fprintf(stderr, "Failed to allocate matrix W.\n");
        free(x->data); // Avoid memory leak
        exit(EXIT_FAILURE);
    }
    w->rows = totalballots;
    w->cols = totalgroups;

    printf("\nInstance Summary:\n");
    printf("Total Votes: %d\n", totalvotes);
    printf("Total Candidates: %d\n", totalcandidates);
    printf("Total Groups: %d\n", totalgroups);
    printf("Total Ballots: %d\n", totalballots);

    printf("\nStarting the randomized instance...\n");
    generateVotes(x, w, &totalvotes, &totalcandidates, &totalgroups, &totalballots);
    printf("The randomized instance is finished!\nWould you like a glimpse of the matrices? (y/n)");

    char choice;
    while (1)
    { // Infinite loop for input validation
        if (scanf(" %c", &choice) == 1)
        {                             // Note the space before %c to skip whitespace
            choice = tolower(choice); // Normalize input to lowercase

            if (choice == 'y')
            {
                printf("\nThe matrix of candidates of dimension (cxb) is:\n");
                printMatrix(x);
                printf("\nThe matrix of demographic groups of dimension (bxg) is:\n");
                printMatrix(w);
                printf("\nPress Enter to continue...\n");
                while (getchar() != '\n')
                    ;      // This would clear word by word the input buffer until it's empty.
                getchar(); // Wait for the user to press a key
                break;
            }
            else if (choice == 'n')
            {
                printf("Alright, proceeding without displaying the matrices.\n");
                break;
            }
            else
            {
                printf("Invalid input. Please enter 'y' or 'n': ");
            }
        }
        else
        {
            printf("Invalid input. Please enter 'y' or 'n': ");
            while (getchar() != '\n')
                ; // This would clear word by word the input buffer until it's empty.
        }
    }
}
