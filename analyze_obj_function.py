import os
import pandas as pd
import numpy as np
import pickle
import time
import random
import matplotlib.pyplot as plt
import seaborn as sns

from group_opt import (
    bootstrap,
    evaluate_all_combs_A,
    group_election_routine,
    new_age_group_name,
    new_group_matrix,
    optimze_age_groups_v2
)

def run_combinations(n_combs=50, output_file='combinations_results.csv', seed=1234):
    """Run the group combination calculations and save results to a CSV."""
    # Set seed for reproducibility
    np.random.seed(seed)
    random.seed(seed)

    # Load data files
    electors = pd.read_csv('2021_11_Presidencial/output/2021_11_Presidencial_ELECTORES.csv')
    votes = pd.read_csv('2021_11_Presidencial/output/2021_11_Presidencial_VOTOS.csv')
    
    with open('2021_11_Presidencial/output/GRUPOS.pickle', 'rb') as handle:
        grupos = pickle.load(handle)
    with open('2021_11_Presidencial/output/CANDIDATOS.pickle', 'rb') as handle:
        candidatos = pickle.load(handle)

    # Merge dataframes on the required keys
    llave_mesa = ['CIRCUNSCRIPCION ELECTORAL', 'MESA']
    df = votes.merge(electors[llave_mesa + grupos], on=llave_mesa, how='inner')

    # Randomly select circunscriptional groups
    circs = np.random.choice(electors['CIRCUNSCRIPCION ELECTORAL'].unique(), size=n_combs, replace=False)
    results_dict = {}

    # Process each circ
    for c in circs:
        print(f'Processing {c}')
        df_c = df.loc[df['CIRCUNSCRIPCION ELECTORAL'] == c]
        x = df_c[candidatos].values
        b = df_c[grupos].values

        for a in range(8, 1, -1):  # Ranges from 8 to 2
            evaluate_all_combs_A(x, b, a, S=50, std_thres=0.05, results_dict=results_dict, circ_name=c)

    # Save the results to a CSV file
    results_df = pd.DataFrame.from_dict(results_dict, orient='index')
    os.makedirs('2021_11_Presidencial/combinations', exist_ok=True)
    results_df.to_csv(os.path.join('combinations',output_file), index=False)
    print(f'Results saved to {output_file}')

def plot_correlations(input_file='combinations_results.csv', obj_func = 'obj_func_prop', output_dir='images'):
    """Generate and save correlation plots based on objective functions."""

    # Load the data
    df = pd.read_csv(input_file)
    df = df[(df['A'] != 1) & (df['A'] != 8)]  # Exclude A=1 and A=8

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

  
    correlations = df.groupby(['A', 'circ']).apply(lambda group: group[obj_func].corr(group['log_likelihood']))
    correlations = correlations.reset_index(name='correlation')

    plt.figure(figsize=(7, 4))
    sns.boxplot(data=correlations, x='A', y='correlation', color='skyblue',
                linecolor='black', linewidth=1.2, width=0.5)
    plt.xlabel('Number of groups ($A$)')
    plt.ylabel(f'Corr(std, log-likelihood)')
    plt.axhline(0, color='grey', linestyle='--')
    output_path = os.path.join(output_dir, f'{obj_func}_corr.pdf')
    plt.savefig(output_path)
    print(f'Saved plot: {output_path}')
    plt.show()

  git config --global user.email "you@example.com"
  git config --global user.name "Your Name"

if __name__ == '__main__':
    run_combinations(n_combs=50)
    plot_correlations(input_file='2021_11_Presidencial/combinations/combinations_results.csv', 
                    obj_func='obj_func_prop', output_dir='images') ## obj function selected (can choose other options)
