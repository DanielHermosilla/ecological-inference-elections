import numpy as np
import pickle
import pandas as pd
import itertools
import matplotlib.pyplot as plt

def load_instance_data(I_list, G_list, M, J, lambda_list, seed_list):
    """Loads data from pickle files for given parameters and computes mean error."""
    instance_dict = {}

    for I in I_list:
        for G in G_list:
            instance_dict[I, G] = []
            for lambda_ in lambda_list:
                lambda_error = []
                wrong_seeds = []
                for s in seed_list:
                    name_of_instance = f"J{J}_M{M}_G{G}_I{I}_lambda{int(lambda_)}"
                    dict_path = f'results/{name_of_instance}/full/{s}.pickle'
                    
                    try:
                        with open(dict_path, 'rb') as f:
                            dict_results = pickle.load(f)
                        
                        p_est = dict_results['p_est']
                        p = dict_results['p']
                        error = np.mean(np.abs(p_est - p))
                        lambda_error.append(error)
                    except FileNotFoundError:
                        wrong_seeds.append(s)
                
                # Log the wrong seeds if any
                if len(wrong_seeds) > 0:
                    wrong_seeds_str = ','.join([str(s) for s in wrong_seeds])
                    print(f'{name_of_instance} {wrong_seeds_str}')

                # Append the mean error for the current lambda
                instance_dict[I, G].append(np.mean(lambda_error))

    return instance_dict

def plot_results(lambda_list, instance_dict, G_list, I_list):
    """Plots the results of the analysis."""
    line_styles = ['-', '--']
    marker_styles = ['^', '^', 'o', 'o']
    fill_styles = [None]
    markerface_colors = [None, 'white']

    for idx, (G, I) in enumerate(itertools.product(G_list, I_list)):
        line_style = line_styles[idx % len(line_styles)]
        marker_style = marker_styles[idx % len(marker_styles)]
        fillstyle = fill_styles[idx % len(fill_styles)]
        markerfacecolor = markerface_colors[idx % len(markerface_colors)]

        plt.plot(lambda_list, instance_dict[(I, G)], label=f'G= {G}, C = {I}', linestyle=line_style,
                 marker=marker_style, markersize=6, fillstyle=fillstyle, markevery=10, markerfacecolor=markerfacecolor)

    plt.legend()
    plt.xlabel('λ')
    plt.ylabel('Mean error')

    # Customize x-axis ticks and labels
    plt.xticks(lambda_list[::5], [f'{l}%' for l in lambda_list[::5]])
    
    # Save and display the plot
    plt.savefig('images/lambda_analysis.pdf', bbox_inches='tight')
    plt.show()

def main():
    """Main function to load data, calculate mean errors, and plot results."""
    # Define parameters
    I_list = [2, 3]  # candidates (C)
    G_list = [2, 3]  # groups (G)
    M = 50           # tables (B)
    J = 100          # people per table (I)
    lambda_list = [5 * i for i in range(0, 21)]
    seed_list = [i for i in range(20)]

    # Load data
    instance_dict = load_instance_data(I_list, G_list, M, J, lambda_list, seed_list)

    # Plot results
    plot_results(lambda_list, instance_dict, G_list, I_list)

if __name__ == "__main__":
    main()
