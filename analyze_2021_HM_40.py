import pandas as pd
import numpy as np

n_mesas = 50

election_name = '2021_11_Presidencial'
candidates = ['EDUARDO ARTES', 'FRANCO ALDO PARISI', 'GABRIEL BORIC', 'JOSE ANTONIO KAST', 'MARCO ENRIQUEZ-OMINAMI', 'NULO BLANCO', 'SEBASTIAN SICHEL', 'YASNA PROVOSTE']
candidates_analyze = ['EDUARDO ARTES', 'FRANCO ALDO PARISI', 'GABRIEL BORIC', 'JOSE ANTONIO KAST', 'MARCO ENRIQUEZ-OMINAMI', 'NULO BLANCO', 'SEBASTIAN SICHEL', 'YASNA PROVOSTE']
# candidates_analyze = ['JOSE ANTONIO KAST']
LABEL_CANDIDATOS = ['G', 'C', 'A', 'B', 'F', 'H', 'D', 'E']


# read 40yo file
output_folder = 'output_40yo' 
df_pais_prob_40yo = pd.read_csv(f'{election_name}/{output_folder}/{election_name}_Probs.csv')
# pivot df_pais_prob_40yo to pass GROUP as columns
df_pais_prob_40yo = df_pais_prob_40yo.pivot(index=['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'MESAS'], columns='GROUP', values=candidates).reset_index()

output_folder = 'output_50yo' 
df_pais_prob_50yo = pd.read_csv(f'{election_name}/{output_folder}/{election_name}_Probs.csv')
df_pais_prob_50yo = df_pais_prob_50yo.pivot(index=['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'MESAS'], columns='GROUP', values=candidates).reset_index()

output_folder = 'output_HM' 
df_pais_prob_HM = pd.read_csv(f'{election_name}/{output_folder}/{election_name}_Probs.csv')
df_pais_prob_HM = df_pais_prob_HM.pivot(index=['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'MESAS'], columns='GROUP', values=candidates).reset_index()

# compute the average difference between the voting probability for JOSE ANTONIO KAST in the 40yo for 18-39 and 40+
# intervalos_mesas = [2, 51, 101, np.Inf]
intervalos_mesas = [n_mesas, np.Inf]
groups_list = [('18-39 HM', '40+ HM'), ('18-49 HM', '50+ HM'), ('H', 'M')]
case_list = ['40yo', '50yo', 'HM']
df_list = [df_pais_prob_40yo, df_pais_prob_50yo, df_pais_prob_HM]
df_diferencias = pd.DataFrame()
df_boxplot = pd.DataFrame()

for k in range(len(groups_list)):
    # k = 0
    group_1, group_2 = groups_list[k]
    df = df_list[k]
    for i in range(len(intervalos_mesas)-1):
        # i = 0
        for c in candidates_analyze:
            # c = 'JOSE ANTONIO KAST'
            mask = (df['MESAS'] >= intervalos_mesas[i]) & (df['MESAS'] < intervalos_mesas[i+1]) & df[(c,'HM')].isna()
            diff = np.mean(df.loc[mask, (c, group_1)] - df.loc[mask, (c, group_2)])
            abs_diff = np.mean(np.abs(df.loc[mask, (c, group_1)] - df.loc[mask, (c, group_2)]))
            df_diferencias_parcial = pd.DataFrame({'CASE': case_list[k], 'MESAS' : [f'{intervalos_mesas[i]}-{intervalos_mesas[i+1]}'], 'CANDIDATO': [c], 'ABS DIFERENCIA': [abs_diff], 'DIFERENCIA': [diff]})
            df_diferencias = pd.concat([df_diferencias, df_diferencias_parcial])
            # df_boxplot_parcial = pd.DataFrame({'CASE': case_list[k], 'MESAS' : [f'{intervalos_mesas[i]}-{intervalos_mesas[i+1]}'], 'CANDIDATO': [c], 'DIFERENCIA': [diff]})
            # df_boxplot = pd.concat([df_boxplot, df_boxplot_parcial])
# df_diferencias.to_csv(f'{election_name}/diferencias_votacion.csv', index=False)
print(df_diferencias)

df_diferencias_2 = df_diferencias.pivot(index=['CANDIDATO'], columns=['CASE', 'MESAS'], values='DIFERENCIA').reset_index()
print("\nDIFERENCIAS")
print(df_diferencias_2)

df_diferencias_3 = df_diferencias.pivot(index=['CANDIDATO'], columns=['CASE', 'MESAS'], values='ABS DIFERENCIA').reset_index()
print("\nDIFERENCIAS ABS")
print(df_diferencias_3)

# hay más diferencia entre 40yo que en 50yo
import seaborn as sns 
import matplotlib.pyplot as plt

for c in candidates:
    df_pais_prob_40yo.drop(columns=[(c, 'HM')], inplace=True)
    df_pais_prob_HM.drop(columns=[(c, 'HM')], inplace=True)

df = df_pais_prob_40yo.merge(df_pais_prob_HM, on=['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'MESAS'])
df.columns
df.shape

for c in candidates_analyze:
    # c = 'JOSE ANTONIO KAST'
    mask = (df['MESAS'] >= intervalos_mesas[i]) & (df['MESAS'] < intervalos_mesas[i+1])
    # df.loc[mask, (c, 'DIF 40')] = df.loc[mask, (c, '40+ HM')] - df.loc[mask, (c, '18-39 HM')]
    # df.loc[mask, (c, 'DIF HM')] = df.loc[mask, (c, 'H')] - df.loc[mask, (c, 'M')]
    df.loc[mask, (c, 'DIF 40')] = np.abs(df.loc[mask, (c, '40+ HM')] - df.loc[mask, (c, '18-39 HM')])
    # df.loc[mask, (c, 'DIF 40')] = np.abs(df.loc[mask, (c, '50+ HM')] - df.loc[mask, (c, '18-49 HM')])
    df.loc[mask, (c, 'DIF HM')] = np.abs(df.loc[mask, (c, 'H')] - df.loc[mask, (c, 'M')])


for c in candidates:
    # df.drop(columns=[(c, '18-39 HM'), (c, '40+ HM'), (c, 'H'), (c, 'M'), (c, 'HM')], inplace=True)
    df.drop(columns=[(c, '18-39 HM'), (c, '40+ HM'), (c, 'H'), (c, 'M')], inplace=True)

mask_row = (df['MESAS'] >= intervalos_mesas[i]) & (df['MESAS'] < intervalos_mesas[i+1]) # & (~df[(c,'DIF 40')].isna())
df = df.loc[mask_row]
df

df.columns = pd.MultiIndex.from_tuples(df.columns)

df = pd.melt(df, 
                  id_vars=[('REGION', ''), ('CIRCUNSCRIPCION ELECTORAL', ''), ('MESAS', '')], 
                  var_name=['Candidate', 'Difference'], 
                  value_name='Value')

# Clean up the column names
df.columns = ['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'MESAS', 'Candidate', 'CASE', 'DIF']

candidate_dict = dict(zip(candidates_analyze, LABEL_CANDIDATOS))
case_dict = {"DIF 40" : "<40, ≥40",
             "DIF HM" : "Men, Women"}

df['CANDIDATOS 2'] = df['Candidate'].replace(candidate_dict)
df['CASE'] = df['CASE'].replace(case_dict)

sns.boxplot(x = df['CANDIDATOS 2'], 
            y = df['DIF'], 
            hue = df['CASE'],
            order=sorted(df['CANDIDATOS 2'].unique()),
            )

plt.xlabel("Candidates", fontsize=12)
plt.ylabel("Absolute Difference of Voting Probability", fontsize=12)
leg = plt.legend(title="Groups:", title_fontsize='10', fontsize='10')  # Customize the legend title and labels
leg._legend_box.align = "left"

plt.savefig(f'images/box_plot_{n_mesas}.pdf', bbox_inches='tight')

# show plot
plt.show()


# /opt/homebrew/bin/python3.9 analize_2021_HM_40.py

# /opt/homebrew/bin/python3.9