#ifndef MEMOIZATION_UTILS_H
#define MEMOIZATION_UTILS_H

#include "uthash.h"
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
    size_t *vector;  // Pointer to the vector `K`
    int vector_size; // Size of the vector
} MemoizationKey;

// Define the memoization table entry
typedef struct
{
    MemoizationKey key;
    double value;
    UT_hash_handle hh; // The hash handle
} MemoizationEntry;

// Memoization table (simple dynamic array for demo purposes)
typedef struct
{
    MemoizationEntry *hashmap;
} MemoizationTable;

// Function prototypes
MemoizationTable *initMemo();
void freeMemo(MemoizationTable *table);
double getMemoValue(MemoizationTable *table, int a, int b, int c, int d, size_t *vector, int vector_size);
void setMemoValue(MemoizationTable *table, int a, int b, int c, int d, size_t *vector, int vector_size, double value);
bool compareKeys(MemoizationKey *key1, MemoizationKey *key2);
uint64_t hashKey(MemoizationKey *key);

#endif // MEMOIZATION_UTILS_H
