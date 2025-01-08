#include "memoizationUtil.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Generates a hash for the key.
 *
 * Given a key structure, it creates the hash for it.
 *
 * @param[in] *key. A pointer towards the key
 *
 * @return A 64-bit hash.
 *
 */
uint64_t hashKey(MemoizationKey *key)
{
    // ---- Initialize the hash as cero
    uint64_t hash = 0;
    for (int i = 0; i < 4; i++)
    { // ---- Associate the hash with indices via a pseudo-random assignation.
        hash ^= key->indices[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    for (int i = 0; i < key->vector_size; i++)
    { // ---- Associate the hash with the elements of the vector used as a key.
        hash ^= (uint64_t)(key->vector[i] * 1e6) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    // ---- Include vector_size for uniqueness, this isn't really necessary for our use case
    hash ^= key->vector_size;
    return hash;
}

/**
 * @brief Initializes a hash table.
 *
 * Building function to initialize a hash table, allocating its memory.
 *
 * @return A pointer towards a MemoizationTable.
 *
 */
MemoizationTable *initMemo()
{
    // ---- Allocates memory for the table
    MemoizationTable *table = (MemoizationTable *)malloc(sizeof(MemoizationTable));
    // ---- Initialize the uthash hash table, initially as NULL.
    table->hashmap = NULL;
    // ---- Return a pointer towards the hash table.
    return table;
}

/**
 * @brief Get the value from the hash table.
 *
 * Given the indices for the hash table, it returns the value stored in the hash table. If there's no value, it'll
 * return "-1.0", whom value is set to "INVALID".
 *
 * @param[in] *table The hash table to query.
 * @param[in] a The first index.
 * @param[in] b The second index.
 * @param[in] c The third index.
 * @param[in] d The fourth index.
 * @param[in] *vector A pointer to the vector used as a key.
 * @param[in] vector_size The size of the vector used as a key.
 *
 * @return double. The value under the key
 *
 */
double getMemoValue(MemoizationTable *table, int a, int b, int c, int d, size_t *vector, int vector_size)
{
    MemoizationEntry *entry;

    // Fill the temporary key structure
    MemoizationKey tempKey;
    tempKey.indices[0] = a;
    tempKey.indices[1] = b;
    tempKey.indices[2] = c;
    tempKey.indices[3] = d;
    tempKey.vector = vector;
    tempKey.vector_size = vector_size;

    // Find using the entire key
    HASH_FIND(hh, table->hashmap, &tempKey, sizeof(MemoizationKey), entry);

    if (entry)
    {
        return entry->value;
    }
    return -1.0; // INVALID
}

/**
 * @brief Create and insert a value into the hash table.
 *
 * Given the indices for the hash table, it inserts a new value of type `double`
 *
 * @param[in, out] *table The hash table to insert the value.
 * @param[in] a The first index.
 * @param[in] b The second index.
 * @param[in] c The third index.
 * @param[in] d The fourth index.
 * @param[in] *vector A pointer to the vector used as a key.
 * @param[in] vector_size The size of the vector used as a key.
 * @param[in] double The value to insert
 *
 * @return void. Changes to be made on the hash table.
 *
 */
void setMemoValue(MemoizationTable *table, int a, int b, int c, int d, size_t *vector, int vector_size, double value)
{
    // ---- Initialize the instance for the new entry
    MemoizationEntry *entry = (MemoizationEntry *)malloc(sizeof(MemoizationEntry));

    // ---- Define each element of the entry structure
    entry->key.indices[0] = a;
    entry->key.indices[1] = b;
    entry->key.indices[2] = c;
    entry->key.indices[3] = d;
    entry->key.vector_size = vector_size;

    // ---- Define the vector used as a key
    entry->key.vector = (size_t *)malloc(vector_size * sizeof(size_t));
    memcpy(entry->key.vector, vector, vector_size * sizeof(size_t));
    // ---- Define the value of the entry
    entry->value = value;

    // ---- Generate the hash key.
    uint64_t keyHash = hashKey(&entry->key);

    // ---- Dinamically add the element with the hash key using HASH_ADD_KEYPTR.
    HASH_ADD_KEYPTR(hh, table->hashmap, &keyHash, sizeof(uint64_t), entry);
}

/**
 * @brief Frees the memory from the Hash Table.
 *
 * Given the hash table, it frees the vector, entry and removes the hash table.
 *
 * @param[in] *table The hash table to be removed
 *
 * @return void
 *
 */
void freeMemo(MemoizationTable *table)
{
    MemoizationEntry *current_entry, *tmp;

    HASH_ITER(hh, table->hashmap, current_entry, tmp)
    {                                            // ---- For each element from the hash table
        free(current_entry->key.vector);         // ---- Free vector memory
        HASH_DEL(table->hashmap, current_entry); // ---- Remove from hash table
        free(current_entry);                     // ---- Free entry memory
    }

    free(table);
}
