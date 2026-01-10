#ifndef UTILS_HASH_EIM
#define UTILS_HASH_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include <stdbool.h>
#include <stdint.h>

#define INVALID -1.0
    // Entry structure: stores the key and the associated value
    typedef struct
    {
        uint64_t key; // Hash of the memoization key
        double value;
    } MemoizationEntry;

    // Hash table structure
    typedef struct
    {
        MemoizationEntry *hashmap; // Pointer to the hash table
    } MemoizationTable;

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
    uint64_t generateHash(int a, int b, int c, int d, const size_t *vector, int vector_size);
    /**
     * @brief Initializes a hash table.
     *
     * Building function to initialize a hash table, allocating its memory.
     *
     * @return A pointer towards a MemoizationTable.
     *
     */
    MemoizationTable *initMemo(void);
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
    double getMemoValue(MemoizationTable *table, int a, int b, int c, int d, size_t *vector, int vector_size);
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
    void setMemoValue(MemoizationTable *table, int a, int b, int c, int d, size_t *vector, int vector_size,
                      double value);
    /**
     * @brief Frees the memory from the Hash Table.
     *
     * Given the hash table, it releases all entries and removes the hash table.
     *
     * @param[in] *table The hash table to be removed
     *
     * @return void
     *
     */
    void freeMemo(MemoizationTable *table);

    /*
     * @brief Creates a hash key for a matrix
     *
     */
    unsigned int computeMatrixKey(const IntMatrix *m);

#ifdef __cplusplus
}
#endif
#endif // HASH_UTILSH
