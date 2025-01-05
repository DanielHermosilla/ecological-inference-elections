#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Define a structure for vector keys
typedef struct
{
    int indices[4];  // The 4 existing indices
    double *vector;  // Pointer to the vector `K`
    int vector_size; // Size of the vector
} MemoizationKey;

// Define the memoization table entry
typedef struct
{
    MemoizationKey key;
    double value;
} MemoizationEntry;

// Memoization table (simple dynamic array for demo purposes)
typedef struct
{
    MemoizationEntry *entries;
    int capacity;
    int size;
} MemoizationTable;

MemoizationTable *initMemo(int capacity)
{
    MemoizationTable *table = (MemoizationTable *)malloc(sizeof(MemoizationTable));
    table->entries = (MemoizationEntry *)calloc(capacity, sizeof(MemoizationEntry));
    table->capacity = capacity;
    table->size = 0;
    return table;
}

// Function to compare two memoization keys
bool compareKeys(MemoizationKey *key1, MemoizationKey *key2)
{
    if (memcmp(key1->indices, key2->indices, 4 * sizeof(int)) != 0)
        return false;
    if (key1->vector_size != key2->vector_size)
        return false;
    return memcmp(key1->vector, key2->vector, key1->vector_size * sizeof(double)) == 0;
}

// Retrieve a value from the table
double getMemoValue(MemoizationTable *table, int a, int b, int c, int d, double *vector, int vector_size)
{
    for (int i = 0; i < table->size; i++)
    {
        // Compare using memcmp
        // 1. If the indices match; and
        // 2. The size of the input vector matches the table vector size; and
        // 3. The vector key is the same as the vector input
        if (memcmp(table->entries[i].key.indices, (int[]){a, b, c, d}, 4 * sizeof(int)) == 0 &&
            table->entries[i].key.vector_size == vector_size &&
            memcmp(table->entries[i].key.vector, vector, vector_size * sizeof(double)) == 0)
        {
            return table->entries[i].value;
        }
    }
    return -1.0; // INVALID
}

// Store a value in the table
void setMemoValue(MemoizationTable *table, int a, int b, int c, int d, double *vector, int vector_size, double value)
{
    if (table->size >= table->capacity)
    {
        table->capacity *= 2;
        table->entries = (MemoizationEntry *)realloc(table->entries, table->capacity * sizeof(MemoizationEntry));
    }
    table->entries[table->size].key.indices[0] = a;
    table->entries[table->size].key.indices[1] = b;
    table->entries[table->size].key.indices[2] = c;
    table->entries[table->size].key.indices[3] = d;
    table->entries[table->size].key.vector_size = vector_size;
    table->entries[table->size].key.vector = (double *)malloc(vector_size * sizeof(double));
    memcpy(table->entries[table->size].key.vector, vector, vector_size * sizeof(double));
    table->entries[table->size].value = value;
    table->size++;
}

// Free the memoization table
void freeMemo(MemoizationTable *table)
{
    for (int i = 0; i < table->size; i++)
    {
        free(table->entries[i].key.vector);
    }
    free(table->entries);
    free(table);
}
