import numpy as np
import pickle
import json
import pandas as pd
import time
import sys
import os

from instance_maker import gen_instance, gen_instance_v2, gen_instance_v3
from EM_full import EM_full, compute_q_list
from EM_mvn_pdf import EM_mvn_pdf, compute_q_mvn_pdf
from EM_simulate import EM_simulate, simulate_Z
from EM_mult import EM_mult, compute_q_multinomial
from EM_mvn_cdf import EM_mvn_cdf, compute_q_mvn_cdf
from EM_algorithm import EM_algorithm, get_p_est

verbose = True
load_bar = True
# convergence_value = 0.001
max_iterations = sys.maxsize
step_size = 1000
python_arg = True

full = lambda X, b, dict_results, dict_file, p_est: EM_full(X, b, p_est = p_est, max_iterations=max_iterations, convergence_value=convergence_value, verbose=verbose, load_bar=load_bar, dict_results = dict_results, save_dict = True, dict_file = dict_file)
EM_methods = [full]
EM_method_names = ["full"]

convergence_values = [0.001, 0.0001]
convergence_names = ['1000', '10000']

J_list = [100] # personas
M_list = [50] # mesas
G_list = [2,3] # grupos
I_list = [2,3] # candidatos 
lambda_list = [0.5]
seed_list = [i+1 for i in range(20)] # seed for starting_point

# start_point_seed = [i+1 for i in range(20)]
# seed_list = [i+1 for i in range(20)]

instances = []

n_instances = len(J_list)*len(M_list)*len(G_list)*len(I_list)*len(seed_list)

for j in J_list:
    for m in M_list:
        for g in G_list:
            for i in I_list:
                for lambda_ in lambda_list:
                    for seed in seed_list:
                        instances.append((j,m,g,i,lambda_, seed))




if __name__ == '__main__':
    if python_arg:
        method_number = int(sys.argv[1])
        instance_number = int(sys.argv[2])
        convergence_value_number = int(sys.argv[3])
        assert (method_number < len(EM_methods)) and (method_number >= 0), f'Method does not exist, should be int between 0 and {len(EM_methods)-1}' 
        assert (instance_number < len(instances)) and (instance_number >= 0), f'Instance does not exist, should be int between 0 and {len(instances)-1}'
        J, M, G, I, lambda_, seed = instances[instance_number]
        EM_method = EM_methods[method_number]
        convergence_value = convergence_values[convergence_value_number]
        method_name = EM_method_names[method_number] + '_cv' + convergence_names[convergence_value_number] + '_convergence'
    else:
        instance_number = np.nan
        J, M, G, I, lambda_, seed = 200, 50, 4, 2, 0.5, 12
        method_number = 3
        method_name = EM_method_names[method_number]
        EM_method = EM_methods[method_number]
    print('-'*70)
    print('instancia ',instance_number,': ',f"J = {J}, M = {M}, G = {G}, I = {I}, lambda = {int(100*lambda_)}%, seed = {seed}")
    print('método ',method_number,': ',EM_method_names[method_number])
    print('-'*70)
    # # generate folder for method if it doesn't exist
    # method_folder = f"results/{method_name}"
    # if not os.path.exists(method_folder):
    #     os.makedirs(method_folder)
    # # generate folder for instance if it doesn't exist
    # instance_folder = f"{method_folder}/J{J}_M{M}_G{G}_I{I}_lambda{lambda_}"
    # if not os.path.exists(instance_folder):
    #     os.makedirs(instance_folder)

    # generate folder for instance if it doesn't exist
    instance_folder = f"results/J{J}_M{M}_G{G}_I{I}_lambda{int(100*lambda_)}"
    if not os.path.exists(instance_folder):
        os.makedirs(instance_folder)
    # generate folder for method_instance if it doesn't exist
    method_instance_folder = f"{instance_folder}/{method_name}"
    if not os.path.exists(method_instance_folder):
        os.makedirs(method_instance_folder)

    # gen instance
    name_of_instance = f"J{J}_M{M}_G{G}_I{I}_L{int(100*lambda_)}"

    gen_instance_v3(G, I, M, J, lambda_ = lambda_, seed = 0, name = name_of_instance, terminar = False)
    # load instance with json  
    with open(f"instances/{name_of_instance}.json", 'r') as f:
        data = json.load(f)
    X = np.array(data["n"]) 
    b = np.array(data["b"])
    p = np.array(data["p"])
    dict_results = {}
    dict_results['X'] = X
    dict_results['b'] = b
    dict_results['p'] = p
    dict_results['J'] = J
    dict_results['M'] = M
    dict_results['G'] = G
    dict_results['I'] = I
    dict_results['lambda_'] = lambda_
    dict_results['seed'] = seed
    dict_results['method'] = method_name
    dict_results['convergence_value'] = convergence_value

    np.random.seed(seed)
    p_est_start = np.random.dirichlet([1]*I, size=G)

    # run EM
    results = EM_method(X, b, dict_results, dict_file = f"{method_instance_folder}/{seed}.pickle", p_est = p_est_start)
    if verbose: print(dict_results['p'])