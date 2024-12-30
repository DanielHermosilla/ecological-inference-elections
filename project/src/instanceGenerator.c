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

int main()
{
    return 1;
}
