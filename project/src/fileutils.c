#include "matrixUtils.h"
#include <stdio.h>
#include <stdlib.h>

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

Matrix *readMatrices(const char *filename, int *count)
{
    FILE *file = fopen(filename, "rb");
    if (!file)
    {
        perror("Failed to open file for reading");
        exit(EXIT_FAILURE);
    }

    // Read the number of matrices
    if (fread(count, sizeof(int), 1, file) != 1)
    {
        perror("Failed to read matrix count");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the matrices
    Matrix *matrices = (Matrix *)malloc(*count * sizeof(Matrix));
    if (!matrices)
    {
        perror("Failed to allocate memory for matrices");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Read each matrix
    for (int i = 0; i < *count; i++)
    {
        if (fread(&matrices[i].rows, sizeof(int), 1, file) != 1 || fread(&matrices[i].cols, sizeof(int), 1, file) != 1)
        {
            perror("Failed to read matrix metadata");
            fclose(file);
            free(matrices);
            exit(EXIT_FAILURE);
        }

        size_t size = matrices[i].rows * matrices[i].cols;
        matrices[i].data = (double *)malloc(size * sizeof(double));
        if (!matrices[i].data)
        {
            perror("Failed to allocate memory for matrix data");
            fclose(file);
            free(matrices);
            exit(EXIT_FAILURE);
        }

        if (fread(matrices[i].data, sizeof(double), size, file) != size)
        {
            perror("Failed to read matrix data");
            fclose(file);
            free(matrices[i].data);
            free(matrices);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
    return matrices;
}
