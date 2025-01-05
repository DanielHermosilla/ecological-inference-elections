#include "uthash.h"
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
    UT_hash_handle hh; // The hash handle
} MemoizationEntry;

// Memoization table (simple dynamic array for demo purposes)
typedef struct
{
    MemoizationEntry *hashmap;
} MemoizationTable;

// Generate a hash for the key
uint64_t hashKey(MemoizationKey *key)
{
    uint64_t hash = 0;
    for (int i = 0; i < 4; i++)
    {
        hash ^= key->indices[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    for (int i = 0; i < key->vector_size; i++)
    {
        hash ^= (uint64_t)(key->vector[i] * 1e6) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    hash ^= key->vector_size; // Include vector_size for uniqueness
    return hash;
}

MemoizationTable *initMemo()
{
    MemoizationTable *table = (MemoizationTable *)malloc(sizeof(MemoizationTable));
    table->hashmap = NULL; // Initialize uthash hash table
    return table;
}

double getMemoValue(MemoizationTable *table, int a, int b, int c, int d, double *vector, int vector_size)
{
    MemoizationEntry *entry;

    MemoizationKey tempKey;
    tempKey.indices[0] = a;
    tempKey.indices[1] = b;
    tempKey.indices[2] = c;
    tempKey.indices[3] = d;
    tempKey.vector = vector;
    tempKey.vector_size = vector_size;

    uint64_t keyHash = hashKey(&tempKey);

    HASH_FIND(hh, table->hashmap, &keyHash, sizeof(uint64_t), entry);

    if (entry)
    {
        return entry->value;
    }
    return -1.0; // INVALID
}

void setMemoValue(MemoizationTable *table, int a, int b, int c, int d, double *vector, int vector_size, double value)
{
    MemoizationEntry *entry = (MemoizationEntry *)malloc(sizeof(MemoizationEntry));

    entry->key.indices[0] = a;
    entry->key.indices[1] = b;
    entry->key.indices[2] = c;
    entry->key.indices[3] = d;
    entry->key.vector_size = vector_size;

    entry->key.vector = (double *)malloc(vector_size * sizeof(double));
    memcpy(entry->key.vector, vector, vector_size * sizeof(double));

    entry->value = value;

    uint64_t keyHash = hashKey(&entry->key);

    HASH_ADD(hh, table->hashmap, keyHash, sizeof(uint64_t), entry);
}

void freeMemo(MemoizationTable *table)
{
    MemoizationEntry *current_entry, *tmp;

    HASH_ITER(hh, table->hashmap, current_entry, tmp)
    {
        free(current_entry->key.vector);         // Free vector memory
        HASH_DEL(table->hashmap, current_entry); // Remove from hash table
        free(current_entry);                     // Free entry memory
    }

    free(table);
}
