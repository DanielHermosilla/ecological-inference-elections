#include <stdio.h>
#include <stdlib.h>

#define INVALID -1 // Invalid value when the element is empty.

typedef struct
{
    double ****data;
    int dim1; // The dimensions
    int dim2;
    int dim3;
    int dim4;
} MemoizationTable;

MemoizationTable *initMemo(int dim1, int dim2, int dim3, int dim4)
{
    MemoizationTable *table = (MemoizationTable *)malloc(sizeof(MemoizationTable));
    table->dim1 = dim1;
    table->dim2 = dim2;
    table->dim3 = dim3;
    table->dim4 = dim4;

    table->data = (double ****)malloc(dim1 * sizeof(double ***));
    for (int i = 0; i < dim1; i++)
    {
        table->data[i] = (double ***)malloc(dim2 * sizeof(double **));
        for (int j = 0; j < dim2; j++)
        {
            table->data[i][j] = (double **)malloc(dim3 * sizeof(double *));
            for (int k = 0; k < dim3; k++)
            {
                table->data[i][j][k] = (double *)malloc(dim4 * sizeof(double));
                for (int l = 0; l < dim4; l++)
                {
                    table->data[i][j][k][l] = INVALID; // Start as an invalid data.
                }
            }
        }
    }
    return table;
}

// Frees all the data from the table.
void freeMemo(MemoizationTable *table)
{
    for (int i = 0; i < table->dim1; i++)
    {
        for (int j = 0; j < table->dim2; j++)
        {
            for (int k = 0; k < table->dim3; k++)
            {
                free(table->data[i][j][k]);
            }
            free(table->data[i][j]);
        }
        free(table->data[i]);
    }
    free(table->data);
    free(table);
}

// Get the values of the table
double getMemoValue(MemoizationTable *table, int a, int b, int c, int d)
{
    if (a < table->dim1 && b < table->dim2 && c < table->dim3 && d < table->dim4)
    {
        return table->data[a][b][c][d];
    }
    return INVALID;
}

// Set the values of the table
void setMemoValue(MemoizationTable *table, int a, int b, int c, int d, int value)
{
    if (a < table->dim1 && b < table->dim2 && c < table->dim3 && d < table->dim4)
    {
        table->data[a][b][c][d] = value;
    }
}
