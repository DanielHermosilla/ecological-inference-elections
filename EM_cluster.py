import numpy as np
import pickle
import json
import pandas as pd
import time
import sys
import os

from instance_maker import gen_instance_v3
from EM_full import EM_full, compute_q_list
from EM_mvn_pdf import EM_mvn_pdf, compute_q_mvn_pdf
from EM_simulate import EM_simulate, simulate_Z
from EM_mult import EM_mult, compute_q_multinomial
from EM_mvn_cdf import EM_mvn_cdf, compute_q_mvn_cdf
from EM_algorithm import EM_algorithm, get_p_est

"""
This script runs Expectation-Maximization (EM) algorithms with user-defined parameters.

Usage:
    python3 EM_cluster.py <method_name> <J> <M> <G> <I> <lambda_> <seed_instance> <seed_pinit>

Arguments:
    - method_name: EM method ('full', 'simulate_100', 'simulate_1000', 'cdf', 'pdf', 'mult').
    - J, M, G, I: Parameters for election instance: number of voters per ballot-box, number of ballot-boxes, number of groups, number of candidates.
    - lambda_: Mixing parameter
    - seed_instance: Seed for generating instances.
    - seed_pinit: Seed for initializing probabilities (-1 for default initialization).

Example:
    python3 <script_name>.py simulate_100 10 5 3 20 0.05 12345 67890

Description:
    Generates instances, runs the chosen EM method, and saves results in the 'results' folder.
"""


verbose = True
load_bar = True
convergence_value = 0.001
max_iterations = sys.maxsize
step_size = 3000
python_arg = True
HR_unique = False

full = lambda X, b, dict_results, dict_file, p_est: EM_full(X, b, p_est=p_est, max_iterations=max_iterations, convergence_value=convergence_value, verbose=verbose, load_bar=load_bar, dict_results = dict_results, save_dict = True, dict_file = dict_file)
simulate_100 = lambda X, b, dict_results, dict_file, p_est: EM_simulate(X, b, p_est=p_est, max_iterations=max_iterations, convergence_value=convergence_value, simulate=True, samples = 100, step_size = step_size, verbose=verbose, load_bar=load_bar, dict_results = dict_results, save_dict = True, dict_file = dict_file, unique = HR_unique)
simulate_1000 = lambda X, b, dict_results, dict_file, p_est: EM_simulate(X, b, p_est=p_est, max_iterations=max_iterations, convergence_value=convergence_value, simulate=True, samples=1000, step_size = step_size, verbose=verbose, load_bar=load_bar, dict_results = dict_results, save_dict = True, dict_file = dict_file, unique = HR_unique)
cdf = lambda X, b, dict_results, dict_file, p_est: EM_mvn_cdf(X, b, p_est=p_est, max_iterations=max_iterations, convergence_value=convergence_value, verbose=verbose, load_bar=load_bar, dict_results = dict_results, save_dict = True, dict_file = dict_file)
pdf = lambda X, b, dict_results, dict_file, p_est: EM_mvn_pdf(X, b, p_est=p_est, max_iterations=max_iterations, convergence_value=convergence_value, verbose=verbose, load_bar=load_bar, dict_results = dict_results, save_dict = True, dict_file = dict_file)
mult = lambda X, b, dict_results, dict_file, p_est: EM_mult(X, b, p_est=p_est, max_iterations=max_iterations, convergence_value=convergence_value, verbose=verbose, load_bar=load_bar, dict_results = dict_results, save_dict = True, dict_file = dict_file)
EM_methods_dic = {'full': full, 'simulate_100': simulate_100, 
                'simulate_1000': simulate_1000,
                'cdf': cdf, 'pdf': pdf, 'mult': mult}


if __name__ == '__main__':
    # print_run_instances()
    # exit()
    if python_arg:
        # read parameters from argv
        method_name = sys.argv[1]
        J = int(sys.argv[2])
        M = int(sys.argv[3])
        G = int(sys.argv[4])
        I = int(sys.argv[5])
        lambda_ = float(sys.argv[6])
        seed_instance = int(sys.argv[7])
        seed_pinit = int(sys.argv[8])

        EM_method = EM_methods_dic[method_name]

        method_name = method_name + '_cv' + str(int(1/convergence_value))
        
    print('-'*70)
    print(f"J = {J}, M = {M}, G = {G}, I = {I}, lambda = {int(100*lambda_)}%, seed_instance = {seed_instance}, seed_pinit = {seed_pinit}")
    print('método ',method_name)
    print('-'*70)


    # generate folder for instance if it doesn't exist
    instance_folder = f"results/J{J}_M{M}_G{G}_I{I}_lambda{int(100*lambda_)}"
    if not os.path.exists(instance_folder):
        os.makedirs(instance_folder)
    # generate folder for method_instance if it doesn't exist
    method_instance_folder = f"{instance_folder}/{method_name}"
    if not os.path.exists(method_instance_folder):
        os.makedirs(method_instance_folder)

    # gen instance
    name_of_instance = f"J{J}_M{M}_G{G}_I{I}_L{int(100*lambda_)}_seed{seed_instance}"

    # create folder instances if it doesn't exist
    if not os.path.exists("instances"):
        os.makedirs("instances")

    gen_instance_v3(G, I, M, J, lambda_ = lambda_, seed = seed_instance, name = name_of_instance, terminar = False)
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
    dict_results['seed'] = seed_instance
    dict_results['method'] = method_name
    dict_results['convergence_value'] = convergence_value

    print(p)

    dict_file_name = f"{method_instance_folder}/{seed_instance}"
    if seed_pinit == -1:
        p_est_start = None
    else:
        # draw from Dirichlet
        np.random.seed(seed_pinit)
        p_est_start = np.random.dirichlet([1]*I, size=G)
        dict_file_name += f"_pinit{seed_pinit}"
    dict_file_name += ".pickle"

    # run EM
    results = EM_method(X, b, dict_results, dict_file = dict_file_name,
                         p_est = p_est_start)

