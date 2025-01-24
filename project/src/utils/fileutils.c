#include "fileutils.h"
#include <cjson/cJSON.h>
#include <errno.h>
#include <libgen.h> // For dirname
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

/**
 * @brief Writes an array of matrices to a binary file, including their metadata.
 *
 * Given a pointer to an array of matrices, it writes them to a binary file. Usually it's useful for saving the `X` and
 * `W` matrices with their probabilities matrix. It's also possible to save a single matrix.
 * @param filename Name of the binary file.
 * @param matrices Pointer to the array of matrices.
 * @param count Number of matrices in the array.
 */
void writeMatrices(const char *filename, Matrix *matrices, int count)
{
    FILE *file = fopen(filename, "wb");
    if (!file)
    {
        perror("Failed to open file for writing");
        exit(EXIT_FAILURE);
    }

    // Write the number of matrices (metadata)
    if (fwrite(&count, sizeof(int), 1, file) != 1)
    {
        perror("Failed to write the number of matrices");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Loop through each matrix and write its metadata and data
    for (int i = 0; i < count; i++)
    {
        Matrix *matrix = &matrices[i];

        // Write matrix metadata (rows and cols)
        if (fwrite(&matrix->rows, sizeof(int), 1, file) != 1 || fwrite(&matrix->cols, sizeof(int), 1, file) != 1)
        {
            perror("Failed to write matrix metadata");
            fclose(file);
            exit(EXIT_FAILURE);
        }

        // Write matrix data
        size_t size = matrix->rows * matrix->cols;
        if (fwrite(matrix->data, sizeof(double), size, file) != size)
        {
            perror("Failed to write matrix data");
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
    printf("Successfully wrote %d matrices to %s\n", count, filename);
}

/**
 * @brief Reads an array of matrices from a binary file.
 *
 * @param filename Name of the binary file.
 * @param matrices Pointer to an array of matrices.
 * @param count Pointer to store the number of matrices read.
 */

void readMatrices(const char *filename, Matrix *matrices, int count)
{
    FILE *file = fopen(filename, "rb");
    if (!file)
    {
        perror("Failed to open file for reading");
        exit(EXIT_FAILURE);
    }
    // Read the number of matrices
    int matrix_count;
    if (fread(&matrix_count, sizeof(int), 1, file) != 1)
    {
        perror("Failed to read matrix count");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Check if the caller-provided array can hold the matrices
    if (count < matrix_count)
    {
        fprintf(stderr, "Provided array can only hold %d matrices, but file contains %d\n", count, matrix_count);
        fclose(file);
        exit(EXIT_FAILURE);
    }
    count = matrix_count; // Update the count with the actual number
    // Read the number of

    // Read each matrix
    for (int i = 0; i < count; i++)
    {
        if (fread(&matrices[i].rows, sizeof(int), 1, file) != 1 || fread(&matrices[i].cols, sizeof(int), 1, file) != 1)
        {
            perror("Failed to read matrix metadata");
            fclose(file);
            freeMatrix(&matrices[i]);
            exit(EXIT_FAILURE);
        }
        // Create matrix with dynamic data allocation
        matrices[i] = createMatrix(matrices[i].rows, matrices[i].cols);

        // Read matrix data
        size_t size = matrices[i].rows * matrices[i].cols;
        if (fread(matrices[i].data, sizeof(double), size, file) != size)
        {
            perror("Failed to read matrix data");
            fclose(file);
            freeMatrix(&matrices[i]); // Free allocated data
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}

void readJSONAndStoreMatrices(const char *filename, Matrix *w, Matrix *x, Matrix *p)
{
    printf("Filename is %s\n", filename);
    // Open the JSON file
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        perror("Unable to open the JSON file");
        exit(EXIT_FAILURE);
    }

    // Get the file size and read its content
    fseek(file, 0, SEEK_END);
    long fileSize = ftell(file);
    rewind(file);

    char *jsonContent = malloc(fileSize + 1);
    if (!jsonContent)
    {
        perror("Memory allocation for JSON content failed");
        fclose(file);
        exit(EXIT_FAILURE);
    }
    fread(jsonContent, 1, fileSize, file);
    jsonContent[fileSize] = '\0'; // Null-terminate the JSON string
    fclose(file);

    // Parse the JSON content
    cJSON *json = cJSON_Parse(jsonContent);
    if (!json)
    {
        fprintf(stderr, "JSON parsing error: %s\n", cJSON_GetErrorPtr());
        free(jsonContent);
        exit(EXIT_FAILURE);
    }

    // Extract arrays "p", "r", "n", and "b", "n" is candidates, "b" is group
    cJSON *pArray = cJSON_GetObjectItem(json, "p");
    cJSON *rOuterArray = cJSON_GetObjectItem(json, "r");
    cJSON *rArray = cJSON_GetArrayItem(rOuterArray, 0); // Move into the first sub-array
    cJSON *nArray = cJSON_GetObjectItem(json, "n");
    cJSON *bArray = cJSON_GetObjectItem(json, "b");

    if (!cJSON_IsArray(pArray) || !cJSON_IsArray(rArray) || !cJSON_IsArray(nArray) || !cJSON_IsArray(bArray))
    {
        fprintf(stderr, "One or more required arrays are missing or invalid in the JSON.\n");
        cJSON_Delete(json);
        free(jsonContent);
        exit(EXIT_FAILURE);
    }

    // Create and populate Matrix structs
    *p = createMatrix(cJSON_GetArraySize(pArray), cJSON_GetArraySize(cJSON_GetArrayItem(pArray, 0)));
    for (int i = 0; i < p->rows; i++)
    {
        cJSON *row = cJSON_GetArrayItem(pArray, i);
        for (int j = 0; j < p->cols; j++)
        {
            MATRIX_AT_PTR(p, i, j) = cJSON_GetArrayItem(row, j)->valuedouble;
        }
    }

    Matrix rMatrix = createMatrix(cJSON_GetArraySize(rArray), cJSON_GetArraySize(cJSON_GetArrayItem(rArray, 0)));
    for (int i = 0; i < rMatrix.rows; i++)
    {
        cJSON *row = cJSON_GetArrayItem(rArray, i);
        for (int j = 0; j < rMatrix.cols; j++)
        {
            MATRIX_AT(rMatrix, i, j) = cJSON_GetArrayItem(row, j)->valuedouble;
        }
    }

    *x = createMatrix(cJSON_GetArraySize(cJSON_GetArrayItem(nArray, 0)), cJSON_GetArraySize(nArray));
    for (int i = 0; i < x->cols; i++)
    {
        cJSON *row = cJSON_GetArrayItem(nArray, i);
        for (int j = 0; j < x->rows; j++)
        {
            MATRIX_AT_PTR(x, j, i) = cJSON_GetArrayItem(row, j)->valuedouble;
        }
    }

    *w = createMatrix(cJSON_GetArraySize(bArray), cJSON_GetArraySize(cJSON_GetArrayItem(bArray, 0)));
    for (int i = 0; i < w->rows; i++)
    {
        cJSON *row = cJSON_GetArrayItem(bArray, i);
        for (int j = 0; j < w->cols; j++)
        {
            MATRIX_AT_PTR(w, i, j) = cJSON_GetArrayItem(row, j)->valuedouble;
        }
    }

    // Print or use matrices as needed
    /*
    printf("Real probabilities matrix:\n");
    printMatrix(p);
    // printf("Matrix r:\n");
    // printMatrix(&rMatrix);
    printf("Candidate matrix:\n");
    printMatrix(x);
    printf("Groups matrix:\n");
    printMatrix(w);
*/
    // Clean up
    // freeMatrix(&p);
    // freeMatrix(&rMatrix);
    // freeMatrix(&x);
    // freeMatrix(&w);
    cJSON_Delete(json);
    free(jsonContent);
}

/**
 * @brief Write the results to a .txt file.
 *
 * This function writes the name of the input JSON file and the contents of the
 * matrix `p` to a specified output text file.
 *
 * @param outputFileName Name of the output .txt file.
 * @param inputFileName Name of the input .json file.
 * @param p The matrix to write to the file.
 */
void writeResults(const char *outputFileName, const char *inputFileName, const char *methodUsed,
                  const double convergence, const double iterations, const double time, const int iterationsMade,
                  const Matrix *pReal, const Matrix *pCalculated, int S, int M, bool hit)
{
    // Ensure the directory for the output file exists
    char *pathCopy = strdup(outputFileName); // Make a copy of the file path
    char *directory = dirname(pathCopy);     // Extract the directory part of the path

    struct stat st = {0};
    if (stat(directory, &st) == -1) // Check if directory exists
    {
        if (mkdir(directory, 0755) != 0 && errno != EEXIST) // Create directory if it doesn't exist
        {
            perror("Unable to create the output directory");
            free(pathCopy);
            exit(EXIT_FAILURE);
        }
    }
    free(pathCopy); // Free the allocated copy

    // Open the output file for writing
    FILE *file = fopen(outputFileName, "w");
    if (!file)
    {
        perror("Unable to open the output file");
        exit(EXIT_FAILURE);
    }

    // Write the input file name
    fprintf(file, "Input JSON File: %s\n", inputFileName);
    if (!hit)
    {
        fprintf(file,
                "Estimated with the %s method and the following parameters:\nConvergence threshold:\t%.7f\nTotal "
                "iterations made:\t%d\nMaximum "
                "amount of "
                "iterations threshold:\t%d\nTime taken:\t%.5f seconds\n",
                methodUsed, convergence, iterationsMade, (int)iterations, time);
    }
    else
    {
        fprintf(file,
                "Estimated with the %s method and the following parameters:\nSamples (S):\t%d\nStep size "
                "(M):\t%d\nConvergence threshold:\t%.7f\nTotal iterations made:\t%d\nMaximum amount of "
                "iterations threshold:\t%d\nTime taken:\t%.5f seconds\n",
                methodUsed, S, M, convergence, iterationsMade, (int)iterations, time);
    }
    // Write the matrix dimensions
    fprintf(file, "\nThe ground truth probability matrix (%dx%d) is:\n", pReal->rows, pReal->cols);

    // Write the matrix data
    for (int i = 0; i < pReal->rows; i++)
    {
        for (int j = 0; j < pReal->cols; j++)
        {
            fprintf(file, "%.6f ", MATRIX_AT_PTR(pReal, i, j)); // Write with 6 decimal precision
        }
        fprintf(file, "\n");
    }

    fprintf(file, "\nThe estimated probability matrix (%dx%d) is:\n", pCalculated->rows, pCalculated->cols);
    // Write the matrix data
    for (int i = 0; i < pCalculated->rows; i++)
    {
        for (int j = 0; j < pCalculated->cols; j++)
        {
            fprintf(file, "%.6f ", MATRIX_AT_PTR(pCalculated, i, j)); // Write with 6 decimal precision
        }
        fprintf(file, "\n");
    }

    // Close the file
    fclose(file);

    printf("Results successfully written to %s\n", outputFileName);
}
