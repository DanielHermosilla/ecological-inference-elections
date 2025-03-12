import itertools
import numpy as np
from scipy.stats import multinomial
from math import comb
import pickle
import os
import struct


# revisar func
def combinations(I_size, b):
    if b == 0:
        return np.zeros(shape=(1, I_size))
    return np.sum(
        list(
            itertools.combinations_with_replacement(np.identity(I_size, dtype=int), b)
        ),
        axis=1,
    )


# para arreglar caso cuando es elegir sobre 0 elementos, revisar
def multinomial_pmf_fixed(x, n, p):
    pmf = multinomial.pmf(x, n, p)
    if np.isnan(pmf):
        return 1.0
    return pmf


def find_tuple(K, k):
    return np.where(np.all(K == k, axis=1))[0][0]


def comb_(x, y):
    if x > 0:
        return comb(x, y)
    else:
        return 0


def row_in_array(myarray, myrow):
    return (myarray == myrow).all(-1).any()


def combinations_filtered(I_size, b, n):
    # print('n',n)
    # print(type(n))
    if b == 0:
        return np.zeros(shape=(1, I_size), dtype=int)
    l_1 = reversed(
        list(
            list(tup)
            for tup in itertools.combinations_with_replacement(range(I_size), b)
        )
    )

    # VERSIÓN ORIGINAL
    l_2 = []
    for tup in l_1:
        tup_2 = [sum([tup[i] == j for i in range(b)]) for j in range(I_size)]
        # for k in range(I_size):
        #     if tup_2[k] > n[k]:
        #         agregar = False
        #         break
        # if agregar:
        #     l_2.append(tup_2)
        if all(np.array(tup_2) <= n):
            l_2.append(tup_2)

    # l_2_ = [[sum([tup[i] == j for i in range(b)]) for j in range(I_size)] for tup in l_1]
    # l_2 = [tup for tup in l_2_ if all(np.array(tup) <= n)]

    # VERSIÓN NUEVA
    # l_2 = np.array([[sum([tup[i] == j for i in range(b)]) for j in range(I_size)]  for tup in l_1])
    # l_2 = l_2[np.all(l_2 <= n,axis=1)]
    # print(l_2)
    # exit()
    return l_2


def vote_transfer(n_m, i_give, i_get, votes):
    n_m[i_give] -= votes
    n_m[i_get] += votes


# write list to binary file
def write_list(a_list, name):

    # store list in binary file so 'wb' mode
    with open(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), "Z_instances", name),
        "wb",
    ) as fp:
        pickle.dump(a_list, fp)
        print("Done writing list into a binary file")


# Read list to memory
def read_list(name):
    # for reading also binary mode is important
    with open(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), "Z_instances", name),
        "rb",
    ) as fp:
        n_list = pickle.load(fp)
        return n_list


def load_election(folder_elections, folder_circunscripcion):
    current_path = os.getcwd()
    election_path = os.path.join(
        current_path, "datos_elecciones", folder_elections, folder_circunscripcion
    )
    b = np.load(os.path.join(election_path, f"b_{folder_circunscripcion}.npy"))
    n = np.load(os.path.join(election_path, f"n_{folder_circunscripcion}.npy"))
    return b, n


# round so that X sums n
def round_sum_n(X):
    n = np.sum(X)
    sorted_index_remainder = np.argsort(X % 1)
    X_round = X // 1
    # print(X_round)
    numbers_to_assign = round(n - np.sum(X_round))
    # print(sorted_index_remainder)
    for i in range(numbers_to_assign):
        X_round[sorted_index_remainder[-(i + 1)]] += 1
    return X_round


def get_partitions(n, A):
    # Get all combinations of A-1 dividers in n-1 gaps
    indices = list(range(1, n))  # Possible indices to place dividers
    partition_list = []
    for dividers in itertools.combinations(indices, A - 1):
        partition = []
        start = 0
        # Add each segment based on dividers
        for divider in dividers:
            partition.append(list(range(start, divider)))
            start = divider
        partition.append(list(range(start, n)))  # Add the last segment
        partition_list.append(partition)
    return partition_list


def readMatrixFromC(filename):
    """
    Reads matrices from a binary file written by the C program.

    Args:
        filename (str): Path to the binary file.

    Returns:
        list of numpy.ndarray: List of NumPy matrices.
    """
    matrices = []

    with open(filename, "rb") as file:
        # Read the number of matrices
        count_data = file.read(4)  # int is 4 bytes in C
        if len(count_data) != 4:
            raise ValueError("Failed to read the number of matrices.")

        count = struct.unpack("i", count_data)[0]  # 'i' -> int
        print(f"Number of matrices: {count}")

        # Read each matrix
        for i in range(count):
            # Read rows and columns (both are ints)
            rows = struct.unpack("i", file.read(4))[0]
            cols = struct.unpack("i", file.read(4))[0]
            print(f"Matrix {i+1}: {rows} x {cols}")

            # Read matrix data
            data_size = rows * cols
            data = file.read(data_size * 8)  # double is 8 bytes in C
            if len(data) != data_size * 8:
                raise ValueError(f"Failed to read data for matrix {i+1}.")

            # Unpack binary data into a NumPy array
            matrix_data = struct.unpack(f"{data_size}d", data)  # 'd' -> double

            # Reshape into a matrix and append to the list
            matrix = np.array(matrix_data).reshape(rows, cols)
            matrices.append(matrix)

    print("All matrices successfully read from binary file.")
    return matrices


if __name__ == "__main__":
    print("Testing get_partitions function")
    # Example usage: Partition the list [0, 1, 2, 3, 4, 5, 6, 7] into A=3 parts
    data_list = list(range(8))  # This is the list [0, 1, 2, 3, 4, 5, 6, 7]
    A = 3
    print(f"All possible ways to partition the list {data_list} into {A} parts:")
    partitions = get_partitions(8, A)
    for i, partition in enumerate(partitions):
        print(f"Partition {i+1}: {partition}")
    print("\n\nTesting readMatrixFromC function")
    matrixList = readMatrixFromC("project/src/matricesTest.bin")
    for i in matrixList:
        print(i)

