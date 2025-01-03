#ifndef MEMOIZATION_UTILS_H
#define MEMOIZATION_UTILS_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INVALID -1.0

// Define a structure for vector keys
typedef struct
{
    int indices[4];  // The 4 existing indices
    double *vector;  // Pointer to the vector (size: TOTAL_CANDIDATES)
    int vector_size; // Size of the vector
} MemoizationKey;

// Define the memoization table entry
typedef struct
{
    MemoizationKey key;
    double value;
} MemoizationEntry;

// Memoization table (hash map-like structure)
typedef struct
{
    MemoizationEntry *entries;
    int capacity;
    int size;
} MemoizationTable;

// Function prototypes
MemoizationTable *initMemo(int capacity);
void freeMemo(MemoizationTable *table);
double getMemoValue(MemoizationTable *table, int a, int b, int c, int d, double *vector, int vector_size);
void setMemoValue(MemoizationTable *table, int a, int b, int c, int d, double *vector, int vector_size, double value);
bool compareKeys(MemoizationKey *key1, MemoizationKey *key2);
uint64_t hashKey(MemoizationKey *key);

#endif // MEMOIZATION_UTILS_H
