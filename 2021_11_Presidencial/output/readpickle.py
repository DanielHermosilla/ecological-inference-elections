import pickle

import pandas as pd

# Load the pickle file
data_dict = pd.read_pickle("2021_11_Presidencial_dfs_distritos.pickle")

# Retrieve the keys from the dictionary
keys_list = list(data_dict.keys())

# Print the keys
print("The keys in the dictionary are:", keys_list)
