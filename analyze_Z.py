import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools

def load_data(G_list, I_list, M, J, step_size, samples, seed_list):
    """Loads data from pickle files for all combinations of G and I."""
    Z_dict = {}
    for G in G_list:
        for I in I_list:
            Z_dict[G, I] = []
            for s in seed_list:
                with open(f"Z_instances/Z_instance_G{G}_I{I}_M{M}_J{J}_step{step_size}_S{samples}_seed{s}_sequence.pickle", 'rb') as f:
                    data = np.array(pickle.load(f))
                    # Flatten the last two dimensions out of four dimensions
                    data = data[:,:,:-1,:-1]
                    data = data.reshape((data.shape[0], data.shape[1], data.shape[2] * data.shape[3]))
                Z_dict[G, I].extend(data)
            Z_dict[G, I] = np.array(Z_dict[G, I])
    return Z_dict

def calculate_mean_correlation(Z_dict, G_list, I_list, S_lim):
    """Calculates the mean correlation for each step size."""
    proms = {}
    for G, I in itertools.product(G_list, I_list):
        Z = Z_dict[G, I]
        proms[G, I] = np.zeros((S_lim, Z.shape[0]))
        for m, Z_s in enumerate(Z):
            for s in range(S_lim):
                proms[G, I][s, m] = np.mean([
                    np.abs(np.corrcoef(Z_s[:-(s+1), i], Z_s[(s+1):, i])[0, 1]) for i in range(Z_s.shape[1])
                ])
        proms[G, I] = np.nanmean(proms[G, I], axis=1)
    return proms

def plot_results(proms, G_list, I_list, S_lim, step_size):
    """Plots the results of the mean correlations."""
    line_styles = ['--', '-']
    marker_styles = ['^', '^', 'o', 'o']
    fill_styles = [None]
    markerface_colors = ['white', None]

    for idx, (G, I) in enumerate(itertools.product(G_list, I_list)):
        line_style = line_styles[idx % len(line_styles)]
        marker_style = marker_styles[idx % len(marker_styles)]
        fillstyle = fill_styles[idx % len(fill_styles)]
        markerfacecolor = markerface_colors[idx % len(markerface_colors)]

        plt.plot(step_size, np.log10(proms[G, I]), label=f'G = {G}, C = {I}', linestyle=line_style,
                 marker=marker_style, markersize=6, fillstyle=fillstyle, markevery=10, markerfacecolor=markerfacecolor)

    plt.legend()
    plt.xlabel('Step size')
    plt.ylabel('Correlation')

    # Customize y-axis ticks and labels
    values_to_show = [1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002]
    y_ticks = [np.log10(v) for v in values_to_show]
    y_ticklabels = [str(v) for v in values_to_show]
    plt.yticks(y_ticks, y_ticklabels)

    # Save and show the figure
    plt.savefig(f'images/corr_graph.pdf', bbox_inches='tight')
    plt.show()

def main():
    """Main function to load data, calculate correlations, and plot results."""
    G_list = [2, 4]
    I_list = [2, 10]
    M = 50
    J = 100
    step_size_value = 100
    samples = 10000
    seed_list = [i for i in range(1, 21)]
    S_lim = 100
    step_size = [i * step_size_value for i in range(1, S_lim + 1)]

    # Load the data
    print('Loading Z instances')
    Z_dict = load_data(G_list, I_list, M, J, step_size_value, samples, seed_list)

    # Calculate mean correlations
    print('Calculating correlations')
    proms = calculate_mean_correlation(Z_dict, G_list, I_list, S_lim)

    # Plot the results
    plot_results(proms, G_list, I_list, S_lim, step_size)
    print('Finished plotting')

if __name__ == "__main__":
    main()