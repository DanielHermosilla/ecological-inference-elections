import numpy as np
import pickle
import matplotlib.pyplot as plt
import os
import pandas as pd

results_folder = "results/"

J_list = [100] # personas
M_list = [50] # mesas
G_list = [2,3] # grupos
I_list = [2,3] # candidatos 
lambda_list = [0.5]
# cv_list = [1000, 10000]
cv_list = [1000]
seed_list = [i+1 for i in range(20)]
# seed_list = [i+1 for i in range(20)]

instances = []
n_instances = len(J_list)*len(M_list)*len(G_list)*len(I_list)*len(seed_list)
EM_method_names = ["full"]

df_differences = []

# read instances
for J in J_list:
    for M in M_list:
        for G in G_list:
            for I in I_list:
                for lambda_ in lambda_list:
                    for cv in cv_list:
                        for seed1 in seed_list:
                            for seed2 in seed_list[seed1:]:
                                pickle_path1 = f'J{J}_M{M}_G{G}_I{I}_lambda{int(100*lambda_)}/full_cv{cv}_convergence/{seed1}.pickle'
                                file_path1 = os.path.join(results_folder, pickle_path1)
                                if os.path.exists(file_path1):
                                    # read pickle
                                    with open(file_path1, 'rb') as handle:
                                        instance1 = pickle.load(handle)
                                    # print(instance)
                                    p_est1 = instance1['p_est']
                                    X1 = instance1['X']
                                else:
                                    print("File not found:", file_path1)
                                # same for seed2
                                pickle_path2 = f'J{J}_M{M}_G{G}_I{I}_lambda{int(100*lambda_)}/full_cv{cv}_convergence/{seed2}.pickle'
                                file_path2 = os.path.join(results_folder, pickle_path2)
                                if os.path.exists(file_path2):
                                    # read pickle
                                    with open(file_path2, 'rb') as handle:
                                        instance2 = pickle.load(handle)
                                    # print(instance)
                                    p_est2 = instance2['p_est']
                                    X2 = instance2['X']
                                # dif_instance = np.max(np.abs(p_est1 - p_est2))
                                dif_instance = np.mean(np.abs(p_est1 - p_est2))
                                # frobenius norm
                                # dif_instance = np.linalg.norm(p_est1 - p_est2)
                                # print(seed1, seed2, dif_instance)
                                # append the dif
                                df_differences.append([J,M,G,I,lambda_,cv,seed1,seed2,dif_instance])
                                # else:
                                #     print("File not found:", file_path)
                                #     print("File path:", file_path)
                                #     print("")

df_differences = pd.DataFrame(df_differences, columns=['J','M','G','I','lambda','cv','seed1','seed2','dif_instance'])


# pasar a latex la tabla
df_differences_latex = df_differences[df_differences['cv'] == 1000]
df_differences_latex = df_differences[['G', 'I', 'dif_instance']]
# round to 4 decimals latex
print(df_differences_latex.groupby(['G','I']).agg(
    ['mean', 'std', 'max'])[['dif_instance']].round(4).to_latex(
        float_format="{:0.4f}".format))