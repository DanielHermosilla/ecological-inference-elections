import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from group_opt import new_group_matrix, optimze_age_groups_v2, bootstrap
from EM_mult import EM_mult
from matplotlib.colors import ListedColormap
import textwrap
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

def load_data(votes_path, electors_path, candidates_path, groups_path):
    votes = pd.read_csv(votes_path)
    electors = pd.read_csv(electors_path)
    with open(candidates_path, 'rb') as f:
        candidatos = pickle.load(f)
    with open(groups_path, 'rb') as f:
        grupos = pickle.load(f)
    return votes, electors, candidatos, grupos

def merge_data(votes, electors, grupos):
    llave_mesa = ['CIRCUNSCRIPCION ELECTORAL', 'MESA']
    return votes.merge(electors[llave_mesa + grupos], on=llave_mesa, how='inner')

def group_heatmap(df_circ, locales, grupos, cmap='Blues', trunc=0.7, output_path=None):
    fig, ax = plt.subplots(2, int(len(locales)/2), figsize=(2*len(locales), 5), sharey=True)
    cbar_ax = fig.add_axes([1.01, 0.3, 0.03, 0.4])

    for i, l in enumerate(locales):
        index = i//int(len(locales)/2), i%int(len(locales)/2)
        group_proportion = df_circ[df_circ['LOCAL'] == l][grupos].T / np.sum(df_circ[df_circ['LOCAL'] == l][grupos].T)
        im = sns.heatmap(group_proportion, ax=ax[index], xticklabels=False, yticklabels=True,
                         cmap=cmap, cbar=True, cbar_ax=cbar_ax, vmin=0, vmax=trunc,
                         cbar_kws={"shrink": 0.3})
        if trunc < 1:
            n_trunc = int(trunc*10)
            cbar_ax.set_yticks([0.1*j for j in range(n_trunc+1)])
            cbar_ax.set_yticklabels([f'   {0.1*j:.1f}' for j in range(n_trunc)] + [f' > {0.1*n_trunc:.1f}'])
        cbar_ax.set_title('Proportion', fontsize=10, pad=10)

        if index[1] == 0:
            ax[index].set_yticklabels(grupos)
            ax[index].set_ylabel('Age Range')
            ax[index].tick_params(left=True)
        else:
            ax[index].tick_params(left=False)

        title = '\n'.join(textwrap.wrap(f'{l}', width=30))
        ax[index].text(0.5, 1.08, title, transform=ax[index].transAxes, va='center', ha='center')

    fig.subplots_adjust(wspace=0.15, right=0.95)

    if output_path:
        fig.savefig(output_path, bbox_inches='tight')
    plt.show()

def plot_probabilities(df, candidatos, grupos, circs, candidates_used, group_combinations_list, output_path=None):
    candidate_label = ['Candidate A', 'Candidate B', 'Candidate C']
    markers = ['o', 'v', 's']
    linestyles = ['-', '--', '-.']

    fig = plt.figure(constrained_layout=True, figsize=(12, 7))
    subfigs = fig.subfigures(nrows=3, ncols=1, height_ratios=[3, 3, 1])

    for i, circ in enumerate(circs):
        x = df.loc[df['CIRCUNSCRIPCION ELECTORAL'] == circ, candidatos].values
        b = df.loc[df['CIRCUNSCRIPCION ELECTORAL'] == circ, grupos].values
        M = len(x)
        subfigs[i].suptitle(f'\n {circ} (B = {M})')
        axs = subfigs[i].subplots(nrows=1, ncols=2, sharey=True)

        for j, group_combinations in enumerate(group_combinations_list):
            b_new, group_names = new_group_matrix(b, group_combinations)
            results = {}
            em = EM_mult(x, b_new, max_iterations=10000, dict_results=results, verbose=False, load_bar=False)
            ll = (results['Q'] - results['E_log_q']) / M
            p_est = results['p_est']

            for c_ix, c in enumerate(candidates_used):
                axs[j].plot(p_est[:, c], label=candidate_label[c_ix], marker=markers[c_ix], linestyle=linestyles[c_ix])

            axs[j].set_xticks(range(len(group_names)), group_names)
            axs[j].text(0.8, 0.9, f'Mean Log-likelihood: {np.round(ll, 3)}',
                        horizontalalignment='center', verticalalignment='center', transform=axs[j].transAxes,
                        bbox=dict(facecolor='white', alpha=0.5))
            axs[0].set_ylim([-0.02, 1.02])
            if j == 0:
                axs[j].set_ylabel('Probability')

    axs_label = subfigs[2].subplots(nrows=1, ncols=1)
    for c_ix, c in enumerate(candidates_used):
        axs_label.plot([0, 0], [0, 0], label=candidate_label[c_ix], marker=markers[c_ix], linestyle=linestyles[c_ix])
    axs_label.set_frame_on(False)
    axs_label.xaxis.set_visible(False)
    axs_label.yaxis.set_visible(False)
    axs_label.legend(loc='center', ncol=3)

    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
    plt.show()

def plot_ll_vs_std(df, candidatos, grupos, circs, group_combinations_list, output_path = None):
    colors = ['green', 'orange', 'purple']
    size_list = [40, 90, 140, 200, 260]

    df_groups = []

    for circ in circs:
        print(circ)
        opt_ll = -np.inf
        x = df.loc[df['CIRCUNSCRIPCION ELECTORAL'] == circ, candidatos].values
        b = df.loc[df['CIRCUNSCRIPCION ELECTORAL'] == circ, grupos].values
        M = len(x)
        # iterate backwards group_cobinations_list
        below_thres = False
        for group_combinations in group_combinations_list[::-1]:
        # group_combinations = optimze_age_groups(b, age_groups = 8, min_std = min_std, max_corr = max_corr)
            b_new, _ = new_group_matrix(b, group_combinations) 
            results = {}
            em = EM_mult(x, b_new, max_iterations = 10000, dict_results = results, verbose = False, load_bar =False)
            ll = results['Q'] - results['E_log_q']
            p_bootstrap = bootstrap(x, b_new, S = 20, seed = 123)
            p_std = np.max(np.std(p_bootstrap, axis = 0))
            if not below_thres:
                if p_std <= 0.05:
                    below_thres = True
                    df_groups.append([circ, M, group_combinations, ll, p_std, 1])
            df_groups.append([circ, M, group_combinations, ll, p_std, 0])    

    df_groups = pd.DataFrame(df_groups, columns = ['Circunscripcion', 'M', 'group_combinations', 'll', 'p_std', 'alphas'])

    # Some formatting for a nice plot
    df_groups['G'] = df_groups['group_combinations'].apply(len) # number of groups
    df_groups['ll-M'] = df_groups['ll']/df_groups['M'] # average ll
    df_groups['Circunscripcion-M'] = df_groups['Circunscripcion'] + ' ' + df_groups['M'].astype(str)
    df_groups['District'] = df_groups['Circunscripcion'].apply(lambda x: x.title())  #.apply(lambda x: x.split(' ')[0])
    df_groups['District'] = df_groups['District'].astype(str) + ' (B = ' + df_groups['M'].astype(str) +')'
    df_groups['# Agg Groups'] = df_groups['G']

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Plot scatter
    scatter = sns.scatterplot(
        data=df_groups, 
        x='ll-M', 
        y='p_std',
        hue='District', 
        size='# Agg Groups', 
        sizes=size_list, 
        legend='brief',  # Enable legend generation
        palette=colors,
        alpha=1, 
        ax=ax
    )

    # Set axis labels
    ax.set_xlabel('Mean Log-likelihood', fontsize=12)
    ax.set_ylabel('Max std(p)', fontsize=12)

    # Set axis limits
    ax.set_ylim([-0.005, 0.55])
    ax.set_xlim([-8, 0.01])

    # Adjust tick parameters
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Get handles and labels for the legend
    h,l = ax.get_legend_handles_labels()
    l1 = ax.legend(h[:len(circs)+1],l[:len(circs)+1], loc='upper left', fontsize = 12)
    l2 = ax.legend(h[len(circs)+1:],l[len(circs)+1:], loc='center left', fontsize = 12)
    ax.add_artist(l1) 

    # Show and save the plot
    plt.show()
    plt.savefig('images/ll_std_example.pdf', bbox_inches='tight')
    plt.close()



def main():
    # File paths
    votes_path = '2021_11_Presidencial/output/2021_11_Presidencial_VOTOS.csv'
    electors_path = '2021_11_Presidencial/output/2021_11_Presidencial_ELECTORES.csv'
    candidates_path = '2021_11_Presidencial/output/CANDIDATOS.pickle'
    groups_path = '2021_11_Presidencial/output/GRUPOS.pickle'

    # Load data
    votes, electors, candidatos, grupos = load_data(votes_path, electors_path, candidates_path, groups_path)
    df = merge_data(votes, electors, grupos)

    # Plot heatmaps
    circ = 'PUENTE ALTO'
    electors_circ = electors[electors['CIRCUNSCRIPCION ELECTORAL'] == circ]
    locales = electors_circ['LOCAL'].unique()
    loc_used = ['COLEGIO PARTICULAR PADRE JOSE KENTENICH', 'COLEGIO NUEVA ERA SIGLO XXI SEDE PUENTE ALTO',
                'LICEO INDUSTRIAL MUNICIPALIZADO A 116 LOCAL: 1', 'COLEGIO MAIPO LOCAL: 2',
                'COLEGIO COMPANIA DE MARIA PUENTE ALTO LOCAL: 2', 'ESCUELA LOS ANDES LOCAL: 1']
    group_heatmap(electors_circ, loc_used, grupos, cmap='Blues', output_path='images/group_distribution_example.pdf')

    # Plot probabilities
    circs = ['EL GOLF', 'SALTOS DEL LAJA']
    group_combinations_list = [[[0], [1], [2], [3], [4], [5], [6], [7]],
                                [[0, 1], [2, 3, 4], [5, 6, 7]]]
    candidates_used = [2, 3, 5]
    plot_probabilities(df, candidatos, grupos, circs, candidates_used, group_combinations_list, 
                       output_path='images/group_probabilities_example.pdf')

    # Plot ll vs std
    circs = ['CANCURA', 'EL BELLOTO', 'PROVIDENCIA']
    group_combinations_list = [[[0,1,2,3,4,5,6,7]],
                    [[0,1,2,3],[4,5,6,7]],
                    [[0,1],[2,3],[4,5],[6,7]],
                    [[0,1],[2],[3],[4],[5],[6,7]],
                    [[0],[1],[2],[3],[4],[5],[6],[7]]]
    plot_ll_vs_std(df, candidatos, grupos, circs, group_combinations_list, 
                   output_path = 'images/ll_vs_std_example.pdf')

if __name__ == '__main__':
    main()




# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# import pickle
# from group_opt import new_group_matrix
# from matplotlib import colors
# import textwrap

# import pandas as pd
# import numpy as np
# from EM_mult import EM_mult
# from group_opt import optimze_age_groups_v2, new_group_matrix, bootstrap
# import pickle
# import warnings
# import matplotlib.pyplot as plt
# import seaborn as sns
# warnings.filterwarnings("ignore", category=RuntimeWarning)


# from matplotlib import colors
# from matplotlib.colors import ListedColormap
# from matplotlib.gridspec import GridSpec

# votes = pd.read_csv('2021_11_Presidencial/output/2021_11_Presidencial_VOTOS.csv')
# electors = pd.read_csv('2021_11_Presidencial/output/2021_11_Presidencial_ELECTORES.csv')
# # read pickles
# with open('2021_11_Presidencial/output/CANDIDATOS.pickle', 'rb') as f:
#     candidatos = pickle.load(f)
# with open('2021_11_Presidencial/output/GRUPOS.pickle', 'rb') as f:
#     grupos = pickle.load(f)

# llave_mesa = ['CIRCUNSCRIPCION ELECTORAL', 'MESA']
# df = votes.merge(electors[llave_mesa + grupos], on = llave_mesa, how = 'inner')

# def group_heatmap(df_circ, locales, cmap = 'Blues', trunc = .7):

#     fig, ax = plt.subplots(2, int(len(locales)/2), figsize=(2*len(locales),5), sharey = True)
                            

#     # max_val = np.max(df_circ[grupos][df_circ['LOCAL'].isin(locales)].T/np.sum(df_circ[grupos][df_circ['LOCAL'].isin(locales)].T))
#     cbar_ax = fig.add_axes([1.01, .3, .03, .4])

#     for i, l in enumerate(locales):
#         index = i//int(len(locales)/2),i%int(len(locales)/2)
#         group_proportion = df_circ[df_circ['LOCAL'] == l][grupos].T/np.sum(df_circ[df_circ['LOCAL'] == l][grupos].T)
#         # group_proportion = group_proportion.values
#         im = sns.heatmap(group_proportion, ax=ax[index], xticklabels=False, yticklabels = True,
#                     cmap = cmap, cbar=True, cbar_ax = cbar_ax ,vmin = 0, vmax = trunc,
#                     cbar_kws={"shrink": 0.3})
#         # colorbar tick label
#         if trunc < 1:
#             n_trunc = int(trunc*10)
#             cbar_ax.set_yticks([0.1*j for j in range(n_trunc+1)])
#             cbar_ax.set_yticklabels([f'   {0.1*j:.1f}' for j in range(n_trunc)] + [f' > {0.1*n_trunc:.1f}'])
#         # cbar title
#         # cbar_ax.set_ylabel('Proportion', rotation=270, labelpad = 15)
#         # cbar title onn top
#         cbar_ax.set_title('Proportion', fontsize = 10, pad = 10)

#         # print((df_circ[df_circ['LOCAL'] == l][grupos].T)/np.sum(df_circ[df_circ['LOCAL'] == l][grupos].T))
#         if index[1] == 0:
#             # ax[index].set_yticks(range(len(grupos)))
#             # ax[index].set_yticklabels(grupos, rotation=90)
#             ax[index].set_yticklabels(grupos)
#             ax[index].set_ylabel('Age Range')
#             ax[index].tick_params(left = True)
#         else:
#             # dont show ticks
#             ax[index].tick_params(left = False)
#         if index[0] == 1:
#             ax[index].set_xlabel('Ballot-box')
#         title_lines = textwrap.wrap(f'{l}', width=30)  # Adjust the width as needed
#         title = '\n'.join(title_lines)
        
#         ax[index].text(0.5, 1.08, title, transform=ax[index].transAxes, va='center', ha='center')
#         ax[index].set_xlabel('')

#     # WSPACE
#     fig.subplots_adjust(wspace = 0.15)
#     # put cbar closer
#     fig.subplots_adjust(right=0.95)

#     # plt.savefig(f'images/elections/group_distribution.pdf')
#     # save fig
#     fig.savefig(f'images/group_distribution_example.pdf', bbox_inches='tight')
#     # # Create the colorbar
#     # cbar = plt.colorbar(im, cax=cbar_ax)

#     # fig.colorbar(im, ax = ax, shrink = 0.8, pad =.01)
#     # plt.tight_layout(w_pad = 0.4)
#     plt.show()


# # seed = 30
# seed  = 119
# print(seed)
# np.random.seed(seed)
# circ = 'PUENTE ALTO'
# electors_circ = electors[electors['CIRCUNSCRIPCION ELECTORAL'] == circ]
# electors_circ = electors_circ[['LOCAL', 'MESA'] + grupos].sort_values(by='LOCAL')
# locales = electors_circ['LOCAL'].unique()
# # loc_used_1 = list(np.sort(np.random.choice(locales, 1, replace = False)))
# # loc_used_2 = list(np.sort(['COLEGIO PARTICULAR PADRE JOSE KENTENICH', 'COLEGIO NUEVA ERA SIGLO XXI SEDE PUENTE ALTO', 'ESCUELA BASICA PADRE ALBERTO HURTADO', 
# # 'LICEO INDUSTRIAL MUNICIPALIZADO A 116 LOCAL: 1', 'COLEGIO MAIPO LOCAL: 2',
# # 'ESCUELA OSCAR BONILLA LOCAL: 1', 
# # 'COLEGIO COMPANIA DE MARIA PUENTE ALTO LOCAL: 2']))
# # loc_used = loc_used_1 + loc_used_2

# loc_used = list(np.sort(['COLEGIO PARTICULAR PADRE JOSE KENTENICH', 'COLEGIO NUEVA ERA SIGLO XXI SEDE PUENTE ALTO', 
# 'LICEO INDUSTRIAL MUNICIPALIZADO A 116 LOCAL: 1', 'COLEGIO MAIPO LOCAL: 2',
# 'COLEGIO COMPANIA DE MARIA PUENTE ALTO LOCAL: 2', 'ESCUELA LOS ANDES LOCAL: 1']))
# group_heatmap(electors_circ, loc_used, cmap = 'Blues')






# # count how many rows each "CIRCUNSCRIPCION ELECTORAL" has
# num_mesas = votes.groupby('CIRCUNSCRIPCION ELECTORAL').count()['MESA']
# # circ_m = num_mesas[(num_mesas >= 20 ) & (num_mesas <= 40)]
# # small_circ = np.random.choice(circ_m.index, size = S_mesas, replace = False)[0]
# # circs = ['EL GOLF', small_circ]
# circs = ['EL GOLF', 'SALTOS DEL LAJA']
# # circs = ['EL GOLF', 'MARCHIGUE']

# group_combinations_list = [[[0],[1],[2],[3],[4],[5],[6],[7]],
#                              [[0,1],[2,3,4],[5,6,7]]]

# prob_df = []

# candidates_used = [2,3,5]

# fig = plt.figure(constrained_layout=True, figsize=(12, 7))
# # fig.suptitle('Figure title')
# subfigs = fig.subfigures(nrows = 3, ncols = 1, height_ratios = [3,3,1])

# # Create subplots

# candidate_label = ['Candidate A', 'Candidate B', 'Candidate C']

# # list of 3 markers
# markers = ['o', 'v', 's']
# # list of 3 linestyles
# linestyles = ['-', '--', '-.']


# for i, circ in enumerate(circs):
#     # opt_ll = -np.inf
#     x = df.loc[df['CIRCUNSCRIPCION ELECTORAL'] == circ, candidatos].values
#     b = df.loc[df['CIRCUNSCRIPCION ELECTORAL'] == circ, grupos].values
#     M = len(x)
#     subfigs[i].suptitle(f'\n {circ} (B = {M})')
#     axs = subfigs[i].subplots(nrows=1, ncols=2, sharey = True)
#     print(circ)
#     for j, group_combinations in enumerate(group_combinations_list):
#     # group_combinations = optimze_age_groups(b, age_groups = 8, min_std = min_std, max_corr = max_corr)
#         b_new, group_names = new_group_matrix(b, group_combinations) 
#         results = {}
#         em = EM_mult(x, b_new, max_iterations = 10000, dict_results = results, verbose = False, load_bar =False)
#         ll = (results['Q'] - results['E_log_q'])/M
#         p_est = results['p_est']
#         prob_df.append([circ, M, group_combinations, ll, p_est])
#         for c_ix, c in enumerate(candidates_used):
#             # plot with markers
#             axs[j].plot(p_est[:,c], label = candidate_label[c_ix], 
#                         marker = markers[c_ix], linestyle = linestyles[c_ix])
#         # set ticks for axis
#         axs[j].set_xticks(range(len(group_names)), group_names)

#         # add ll to plot
#         # axs[j].text(0.5, 0.9, f'Mean Log-likelihood: {np.round(ll, 3)}',
#         #  horizontalalignment='center', verticalalignment='center', transform=axs[j].transAxes)

#         # add ll to plot
#         axs[j].text(0.8, 0.9, f'Mean Log-likelihood: {np.round(ll, 3)}',
#          horizontalalignment='center', verticalalignment='center', transform=axs[j].transAxes,
#          bbox=dict(facecolor='white', alpha=0.5))

#         axs[0].set_ylim([-0.02,1.02])
#         # name of y axis
#         if j == 0:
#             axs[j].set_ylabel('Probability')

#         # if i == 1 and j == 0:
#         #     axs[j].legend(loc = (0.3,-0.8))
# axs_label = subfigs[2].subplots(nrows=1, ncols=1)
# for c_ix, c in enumerate(candidates_used):
#     axs_label.plot([0,0],[0,0], label = candidate_label[c_ix], 
#                    marker = markers[c_ix], linestyle = linestyles[c_ix])
# axs_label.set_frame_on(False)
# axs_label.xaxis.set_visible(False)
# axs_label.yaxis.set_visible(False)
# axs_label.legend(loc = 'center', ncol = 3)
# legend = axs_label.legend(loc = 'center', ncol = 3) 
# legend.get_frame().set_alpha(1)
# # save figure

# # set hspace
# # fig.subplots_adjust(hspace=0.3, wspace = 0.1)

# plt.show()
# plt.savefig('images/group_probabilities_example.pdf', bbox_inches='tight')







# import seaborn as sns
# import matplotlib.pyplot as plt


