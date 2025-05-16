/*
Copyright (c) 2025 fastei team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "utils_hash.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

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

uint64_t generateHash(int a, int b, int c, int d, const size_t *vector, int vector_size)
{
    uint64_t hash = 0xcbf29ce484222325; // FNV offset basis

    // Hash the 4 indices
    hash ^= (uint64_t)a;
    hash *= 0x100000001b3;
    hash ^= (uint64_t)b;
    hash *= 0x100000001b3;
    hash ^= (uint64_t)c;
    hash *= 0x100000001b3;
    hash ^= (uint64_t)d;
    hash *= 0x100000001b3;

    // Hash the vector contents
    for (int i = 0; i < vector_size; i++)
    {
        hash ^= (uint64_t)vector[i];
        hash *= 0x100000001b3;
    }

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
MemoizationTable *initMemo(void)
{
    // ---- Allocates memory for the table
    MemoizationTable *table = (MemoizationTable *)Calloc(1, MemoizationTable);
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
    uint64_t hash = generateHash(a, b, c, d, vector, vector_size);

    // Find the entry
    MemoizationEntry *entry;
    HASH_FIND(hh, table->hashmap, &hash, sizeof(uint64_t), entry);

    if (entry)
    {
        return entry->value;
    }

    return -1.0; // Not found
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
    uint64_t hash = generateHash(a, b, c, d, vector, vector_size);

    // Check if the hash already exists
    MemoizationEntry *entry;
    HASH_FIND(hh, table->hashmap, &hash, sizeof(uint64_t), entry);

    if (entry)
    {
        // Update the existing value
        entry->value = value;
        return;
    }

    // Create a new entry
    entry = Calloc(1, MemoizationEntry);
    entry->hash = hash;
    entry->value = value;

    // Add to the hash table
    HASH_ADD(hh, table->hashmap, hash, sizeof(uint64_t), entry);
}

void deleteEntry(MemoizationTable *table, int a, int b, int c, int d, size_t *vector, int vector_size)
{
    uint64_t hash = generateHash(a, b, c, d, vector, vector_size);

    // Find and delete the entry
    MemoizationEntry *entry;
    HASH_FIND(hh, table->hashmap, &hash, sizeof(uint64_t), entry);

    if (entry)
    {
        HASH_DEL(table->hashmap, entry);
        Free(entry); // Free the entry memory
    }
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
    MemoizationEntry *entry, *tmp;
    HASH_ITER(hh, table->hashmap, entry, tmp)
    {
        HASH_DEL(table->hashmap, entry); // Remove from the hash table
        Free(entry);                     // Free the entry
    }
    Free(table); // Free the table structure itself
}

/*
 * @brief Creates a hash key for a matrix
 *
 */
unsigned int computeMatrixKey(const Matrix *m)
{
    unsigned int hash = 5381;
    hash = ((hash << 5) + hash) ^ (unsigned int)m->rows;
    hash = ((hash << 5) + hash) ^ (unsigned int)m->cols;

    int total = m->rows * m->cols;
    int *data = (int *)m->data;

    // Constantes primas para mezclar posiciones (estas pueden ajustarse)
    const unsigned int prime_row = 73856093;
    const unsigned int prime_col = 19349663;

    for (int i = 0; i < total; i++)
    {
        // Se obtiene la fila y columna a partir del índice
        int row = i / m->cols;
        int col = i % m->cols;
        // Combina el valor, la fila y la columna en una única mezcla
        unsigned int value = (unsigned int)data[i] ^ (row * prime_row) ^ (col * prime_col);
        hash = ((hash << 5) + hash) ^ value;
    }
    return hash;
}
