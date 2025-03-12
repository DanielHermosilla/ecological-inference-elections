#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.stats import binom
import json
from group_opt import new_group_matrix, bootstrap

def load_data(election_name, output_folder, threshold):
    """Load primary election data."""
    df_pais = pd.read_csv(f'{election_name}/{output_folder}/{election_name}_PAIS.csv')
    df_lowp = df_pais[df_pais['LOG P-VALOR'] <= threshold]
    return df_pais, df_lowp[['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'LOCAL']].drop_duplicates().to_numpy()

def plot_ballotbox_histogram(df_pais, image_path, truncate = 500):
    """Plot histogram of ballot boxes by district and save as PDF."""

    df_ballotbox = df_pais.copy()
    df_ballotbox = df_ballotbox.drop_duplicates(subset=['CIRCUNSCRIPCION ELECTORAL'])
    df_ballotbox = df_ballotbox.groupby('NUM MESAS').count()[['LOCAL']]
    df_ballotbox = df_ballotbox.reset_index()
    df_ballotbox.rename(columns={'LOCAL': 'COUNT'}, inplace=True)  # Rename the 'LOCAL' column to 'COUNT'

    if truncate is not None:
        # truncate for visualization
        df_ballotbox.loc[df_ballotbox['NUM MESAS'] > truncate, 'NUM MESAS'] = truncate + 10 

        
    # Plot histogram
    bins = np.arange(0, truncate + 30, 25)
    plt.figure(figsize=(8, 5))
    
    # Use 'COUNT' for the heights of histogram bars
    plt.hist(df_ballotbox['NUM MESAS'], bins=bins, weights=df_ballotbox['COUNT'], 
            edgecolor='k', alpha=0.7, color = 'skyblue')
    
    # Set axis labels and customize x-ticks
    plt.xlabel('Number of Ballot-Boxes')
    plt.ylabel('Number of Districts')
    plt.xticks(np.arange(0, truncate + 100, 100).tolist(), [str(i) for i in np.arange(0, truncate, 100).tolist()] + [f'{truncate}+'])

    # Save the plot to the specified output file
    # Disable grid and remove spines
    plt.grid(False)  # Turn off the grid

    plt.savefig(image_path, bbox_inches='tight')
    plt.show()
    plt.close()





def calculate_pvalue_ranges(df_pais, image_path, log_pvalue_precision, write = False):
    """Calculate p-value ranges and export LaTeX table."""
    range_names = ['$(10^{-\\infty}\\;$&,&$\\;10^{' + str(log_pvalue_precision + 1) + '}]$']
    range_names += ['$(10^{' + str(i) + '}\\;$&,&$\\;10^{' + str(i + 1) + '}]$' for i in range(log_pvalue_precision + 1, 0)]
    df_pais['RANGE'] = pd.cut(df_pais['LOG P-VALOR'], bins=[-np.inf] + [i for i in range(log_pvalue_precision + 1, 1)], labels=range_names, include_lowest=True)

    resumen_pvalores = df_pais.groupby(['RANGE']).count()[['LOCAL']]
    resumen_pvalores = resumen_pvalores.iloc[::-1]
    resumen_pvalores['%'] = (resumen_pvalores['LOCAL'] / resumen_pvalores['LOCAL'].sum() * 100).round(2)
    resumen_pvalores.reset_index(inplace=True)
    
    latex_str = resumen_pvalores.to_latex(index=False, header=['p-value Range', 'Total', '\%'], column_format='ccc', float_format="%.2f")
    
    print(latex_str)

    if write:
        with open(image_path, 'w') as f:
            f.write(latex_str)


def calculate_probabilities(df_pais, circs, circ_labels, election_name, 
                            output_folder, image_path, anonymous=True):
    """Calculate and save voting probabilities heatmap with optional candidate anonymization."""
    
    # Load candidate names
    with open(f'{election_name}/{output_folder}/CANDIDATOS.pickle', 'rb') as handle:
        CANDIDATOS = np.array(pickle.load(handle))
    
    # Define custom display order for candidates
    order_candidatos = [2, 3, 1, 6, 7, 4, 0, 5]  # Adjust this list based on required display order
    
    # Define labels based on anonymity setting
    if anonymous:
        LABEL_CANDIDATOS = [chr(65 + i) for i in range(len(CANDIDATOS))]  # A, B, C, ...
    else:
        LABEL_CANDIDATOS = CANDIDATOS[order_candidatos]

    # Adjust specific candidate names if not anonymized
    if not anonymous:
        LABEL_CANDIDATOS = np.where(LABEL_CANDIDATOS == 'NULO BLANCO', 'NULL/BLANK', LABEL_CANDIDATOS)
        LABEL_CANDIDATOS = np.where(LABEL_CANDIDATOS == 'MARCO ENRIQUEZ-OMINAMI', 'MARCO E.-OMINAMI', LABEL_CANDIDATOS)
        LABEL_CANDIDATOS = np.where(LABEL_CANDIDATOS == 'JOSE ANTONIO KAST', 'JOSÉ A. KAST', LABEL_CANDIDATOS)
        LABEL_CANDIDATOS = np.where(LABEL_CANDIDATOS == 'FRANCO ALDO PARISI', 'FRANCO PARISI', LABEL_CANDIDATOS)

    # Gather voting probability data for specified districts
    prob_list, mesas = [], []
    for circ in circs:
        region = df_pais[df_pais['CIRCUNSCRIPCION ELECTORAL'] == circ]['REGION'].iloc[0]
        mesas.append(df_pais[df_pais['CIRCUNSCRIPCION ELECTORAL'] == circ]['NUM MESAS'].iloc[0])
        probs = pd.read_csv(f'{election_name}/{output_folder}/{region}/{circ}/P_{circ}.csv', index_col=0).T.round(2)
        prob_list.append(probs)

    # Set up the plot
    width_ratios = [prob.shape[1] for prob in prob_list]
    fig, ax = plt.subplots(1, len(circs), figsize=(18, 4), width_ratios=width_ratios)
    
    for ix, df_probs in enumerate(prob_list):
        # Adjust candidate order and anonymized labels
        df_probs = df_probs.iloc[order_candidatos]
        
        ax[ix].set_title(f'{circ_labels[ix]} ($B={mesas[ix]}$)', fontsize=14)
        sns_heatmap = sns.heatmap(df_probs, ax=ax[ix], cmap='YlOrRd', cbar=False, vmin=0, vmax=0.6,
                                  annot=True, annot_kws={"size": 14}, fmt='.2f')
        
        # Set y-axis labels with candidate names (anonymized or real)
        if ix == 0:
            ax[ix].set_yticklabels(LABEL_CANDIDATOS, rotation=0, fontsize=14)
            ax[ix].set_ylabel('Candidates', fontsize=14, labelpad=10)
        else:
            ax[ix].set_yticklabels([])
            
        ax[ix].set_xlabel('Age Ranges', fontsize=14, labelpad=10)
        
        if ix > 0:
            ax[ix].tick_params(axis='y', which='both', length=0)
    
    # Add shared colorbar
    cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])
    colorbar = fig.colorbar(sns_heatmap.collections[0], cax=cbar_ax)
    colorbar.set_label('Voting probability', rotation=-90, fontsize=14, labelpad=30)

    # Save and show the plot
    plt.grid(False)
    plt.gca().set_facecolor('white')
    plt.savefig(image_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)




def plot_group_stats(df_pais, election_name, output_folder, image_path):
    """Calculate aggregated group statistics and plot histogram without saving intermediate files."""

    sns.set(style="whitegrid")

    # Extract unique regions and electoral districts
    circs = df_pais[['REGION', 'CIRCUNSCRIPCION ELECTORAL']].drop_duplicates().to_numpy()
    df_circs = df_pais.groupby(['REGION', 'CIRCUNSCRIPCION ELECTORAL']).count()[['MESA']].reset_index()
    n_circs = len(circs)
    
    # Define age groups and initialize count array
    age_groups = ['18-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+']
    age_group_count = np.zeros(len(age_groups), dtype=int)

    # Dictionary to store counts of groups
    group_count_dict = {}

    # Process each region and electoral district
    for region, circ in circs:
        try:
            # Load group data for each district
            with open(f'{election_name}/{output_folder}/{region}/{circ}/GRUPOS.pickle', 'rb') as handle:
                GRUPOS_CIRC = pickle.load(handle)
            
            # Clean group names and count occurrences
            stripped_groups = [g.strip() for g in GRUPOS_CIRC]
            for ix, age_group in enumerate(age_groups):
                if age_group in stripped_groups:
                    age_group_count[ix] += 1
            for group in stripped_groups:
                group_count_dict[group] = group_count_dict.get(group, 0) + 1
            
            # Update the number of groups for each district in df_circs
            df_circs.loc[(df_circs['REGION'] == region) & 
                        (df_circs['CIRCUNSCRIPCION ELECTORAL'] == circ), 'GRUPOS'] = len(GRUPOS_CIRC)
        except:
            print(f'{region} {circ} NOT FOUND') # This can happen with name mismatches
            n_circs -= 1 # don't count mismatches

    # Convert group counts to DataFrame for plotting
    age_group_proportion = age_group_count / n_circs

    # plot barplot, age_group_proportion vs age_groups
    fig, ax = plt.subplots(1, 1, figsize=(7,5))
    ax.bar(age_groups, age_group_proportion,  
            edgecolor='k', alpha = 0.7, color = 'skyblue')

    # axis name
    ax.set_ylabel('Fraction of Districts')
    ax.set_xlabel('Age Range')


    plt.grid(False)
    plt.gca().set_facecolor('white')
    plt.savefig(os.path.join(image_path, 'group_stats_histograms.pdf'), bbox_inches='tight')
    plt.show()
    plt.close()


    ### Plot boxplot
    fig, axs = plt.subplots(figsize = (7,5))
    sns.boxplot(data = df_circs[['GRUPOS', 'MESA']].dropna().astype(int), 
                x = 'GRUPOS', y = 'MESA', color = 'skyblue',
                ax = axs, linewidth=1.2, linecolor='black')
    # namer of axis
    axs.set_ylabel('Number of ballot-boxes (B)')
    axs.set_xlabel('Number of Aggregated Groups (A)')
    plt.grid(False)
    plt.gca().set_facecolor('white')
    # save fig in images/elections directory
    plt.savefig(os.path.join(image_path,'boxplot_groups.pdf'), bbox_inches='tight')
    plt.show()
    plt.close()


    # Plot barplot
    group_number = df_circs.groupby(['GRUPOS']).count()[['MESA']]
    group_number = group_number.reset_index()
    group_number['proportion'] = np.round(group_number['MESA'] / group_number['MESA'].sum(), 2)
    group_number[['GRUPOS']] = group_number[['GRUPOS']].astype(int) 

    # plot barplot, age_group_proportion vs age_groups
    fig, ax = plt.subplots(1, 1, figsize = (7,5))
    sns.barplot(data = group_number, x = 'GRUPOS', y = 'proportion', 
                edgecolor='black', linewidth=1.2, color = 'lightskyblue')
    # axis name
    ax.set_ylabel('Proportion of Districts')
    ax.set_xlabel('Number of Aggregated Groups (A)')
    plt.grid(False)
    plt.gca().set_facecolor('white')
    plt.savefig(os.path.join(image_path,'barplot_groups.pdf'), bbox_inches='tight')
    plt.show()
    plt.close()



    with open(f'{election_name}/{output_folder}/group_count_dict.json', 'w') as fp:
        json.dump(group_count_dict, fp)





def plot_district_heatmap(df_pais, low_p, election_name, output_folder, image_dir, max_pval=8, seed=32):
    """Generate and save heatmap for selected districts based on voting data and p-values."""


    # Save current rcParams
    original_rc_params = plt.rcParams.copy()
    # Set temporary font size
    plt.rcParams.update({'font.size': 8})

    image_dir = os.path.join(image_dir, 'heatmaps')
    os.makedirs(image_dir, exist_ok = True)

    order_candidatos = [2, 3, 1, 6, 7, 4, 0, 5]  # Adjust this list based on required display order


    # Loop over specified regions and districts
    for region, circ, _ in low_p:
        # print('Plotting: ', circ)
        try:
            # Load candidate and group data
            with open(f'{election_name}/{output_folder}/CANDIDATOS.pickle', 'rb') as handle:
                CANDIDATOS = np.array(pickle.load(handle))
                CANDIDATOS = CANDIDATOS[order_candidatos]
            with open(f'{election_name}/{output_folder}/GRUPOS.pickle', 'rb') as handle:
                GRUPOS = pickle.load(handle)
            
            # Define LABELS for anonimity
            LABEL_CANDIDATOS = [chr(65 + i) for i in range(len(CANDIDATOS))]
        


            # Color map for differences
            map_dif = sns.diverging_palette(-50, 130, s=500, l=75, sep=25, center='light', as_cmap=True)

            # Load p-values and district data
            df_region = pd.read_csv(f'{election_name}/{output_folder}/{region}/{region}.csv')
            df_region_pval = df_region[['MESA', 'CIRCUNSCRIPCION ELECTORAL', 'P-VALOR']].copy()
            df_circ = pd.read_csv(f'{election_name}/output/{region}/{circ}/{circ}.csv')
            prob_circ = pd.read_csv(f'{election_name}/{output_folder}/{region}/{circ}/P_{circ}.csv')
            
            with open(f'{election_name}/output/{region}/{circ}/GRUPOS.pickle', 'rb') as handle:
                grupos = pickle.load(handle)
            
            # Find smallest p-value in the district
            df_circ_pval = df_region_pval[df_region_pval['CIRCUNSCRIPCION ELECTORAL'] == circ]
            mesa_smallest = df_circ_pval[df_circ_pval['P-VALOR'] == df_circ_pval['P-VALOR'].min()].iloc[0]
            
            # Sample random subset and add mesa with smallest p-value if missing
            df_heatmap = df_circ.sample(10, random_state=seed)
            if mesa_smallest['MESA'] not in df_heatmap['MESA'].values:
                df_heatmap = pd.concat([df_heatmap, df_circ[df_circ['MESA'] == mesa_smallest['MESA']]], ignore_index=True)
            
            # Prepare data for heatmaps
            df_heatmap = df_heatmap.drop(columns=['P-VALOR'])
            df_heatmap = df_heatmap.merge(df_region_pval, on=['MESA', 'CIRCUNSCRIPCION ELECTORAL'], how='inner')
            dif_CANDIDATOS = np.array([f'DIF_{c}' for c in CANDIDATOS])
            df_heatmap.loc[df_heatmap['P-VALOR'] == 0, 'P-VALOR'] = 10 ** (max_pval-1)
            # zero_text = f'<{max_pval}'

            # Set up plot dimensions
            w_1, w_2 = 0.4, 0.3
            fig, ax = plt.subplots(1, 4, sharey=True, figsize=(w_1 * (2 * len(CANDIDATOS) + len(grupos) + 2), w_2 * (2 + len(df_heatmap))), 
                                width_ratios=[len(CANDIDATOS), len(CANDIDATOS), len(grupos), 2])

            # Title and layout adjustments
            # fig.suptitle(f'REGION {region} - {circ}', fontsize=8) won't include title for paper
            plt.subplots_adjust(top=len(df_heatmap)/(2 + len(df_heatmap)), wspace=0.1)

            # Index by MESA
            df_heatmap = df_heatmap.set_index('MESA')
            
            # Heatmap 1: Votes
            ax[0].set_title('Votes', fontsize=8)
            sns.heatmap(df_heatmap[CANDIDATOS].astype(int), cmap='Reds', cbar=False, fmt='g',
                        annot=True, annot_kws={'fontsize': 8}, ax=ax[0],
                        xticklabels=[l for l in LABEL_CANDIDATOS])
            ax[0].tick_params(axis='y', rotation=0, labelsize=8)
            ax[0].tick_params(axis='x', labelsize = 8)
            ax[0].set_ylabel('Ballot-box id', fontsize = 8)
            ax[0].set_xlabel('Candidates', fontsize=8)

            # Heatmap 2: Votes minus Expected Votes
            ax[1].set_title('Votes minus Expected Votes', fontsize=8)
            sns.heatmap(df_heatmap[dif_CANDIDATOS].astype(int), cmap=map_dif, cbar=False, fmt='g',
                        annot=True, annot_kws={'fontsize': 8}, ax=ax[1], center=0,
                        xticklabels=[l for l in LABEL_CANDIDATOS])
            ax[1].set_xlabel('Candidates', fontsize=8)
            ax[1].tick_params(axis='x', labelsize = 8)


            # Heatmap 3: Voters by Age Group
            ax[2].set_title('Voters', fontsize=8)
            sns.heatmap(df_heatmap[grupos].astype(int), cmap='Blues', cbar=False, fmt='g',
                        annot=True, annot_kws={'fontsize': 8}, ax=ax[2])
            ax[2].set_xlabel('Age Ranges', fontsize=8)
            ax[2].tick_params(axis='x', labelsize = 8)


            # Heatmap 4: p-value
            ax[3].set_title('p-value', fontsize=8)
            sns.heatmap(df_heatmap[['P-VALOR']], cmap='Purples_r', cbar=False,
                        annot=True, fmt='.1e', annot_kws={'fontsize': 8}, ax=ax[3],
                        xticklabels=False, norm=colors.LogNorm(vmin=10**(max_pval - 1), vmax=1))
            
            # Update text to show logarithm values
            for text in ax[3].texts:
                # Get the original text value and convert to a float
                value = float(text.get_text())
                # Convert to base-10 logarithm and set new text
                log_value = np.log10(value)
                text.set_text(f'{log_value:.2f}') 
                if log_value == max_pval - 1:
                    text.set_text('<'+f'{log_value+1:.2f}')



            # Remove yticks for non-leftmost axes
            for ix in range(1, 4):
                ax[ix].tick_params(left=False)
                ax[ix].set_ylabel('')
            
            # Save and close figure
            plt.savefig(os.path.join(image_dir, f'{circ}-{seed}.pdf'), bbox_inches='tight')
            # plt.show()
            plt.close(fig)
        
        except:
            print(f'{region} {circ} NOT FOUND') # mismatch cases
    
    # Reset rcParams to original values
    plt.rcParams.update(original_rc_params)

# Main function
def main():

    # Global settings
    election_name = '2021_11_Presidencial'
    output_folder = 'output'
    images_dir = 'images'
    threshold = -8
    log_pvalue_precision = -9

    # create directories
    os.makedirs(images_dir, exist_ok = True)
    os.makedirs(os.path.join(images_dir, "elections"), exist_ok = True)

    # Load data
    df_pais, df_lowp = load_data(election_name, output_folder, threshold)

    # Perform analyses
    plot_ballotbox_histogram(df_pais, image_path=os.path.join(images_dir, 'ballotboxes.pdf'))
    calculate_pvalue_ranges(df_pais, image_path=os.path.join(images_dir, 'pvalue_ranges.tex'),
                            log_pvalue_precision = log_pvalue_precision, write = False)
    calculate_probabilities(df_pais, circs=['PROVIDENCIA', 'PLAZA EGANA', 'VALLENAR', 'FRUTILLAR'],
                            circ_labels=['Providencia', 'Plaza Egaña', 'Vallenar', 'Frutillar'],
                            election_name=election_name, output_folder=output_folder,
                            image_path=os.path.join(images_dir, 'probability_results.pdf'))
    plot_group_stats(df_pais, election_name, output_folder, images_dir)
    plot_district_heatmap(df_pais, df_lowp, election_name, output_folder, images_dir, 
                            max_pval=threshold, seed=123)

if __name__ == "__main__":
    main()
