# codigo que analiza la heterogeneidad de los grupos demográficos en las mesas electorales
# usando el test chi-cuadrado sobre tabla de contingencia entre grupos demográficos y mesas
# crea un gráfico con histograma acumulado de distritos según p-valores

import pandas as pd
import time
import numpy as np
from unidecode import unidecode
import re
from scipy.stats import chi2_contingency # test chi-2 con matriz de contingencia
import matplotlib.pyplot as plt



verbose = True # print outputs

# function that preprocess the data of number of voters from each demographical group 
# for each ballot box
def count_groups_per_ballot(election_name):
    # election_name = "2021_11_Presidencial"
    # function that correct names
    def correct_names(col_name):
        if col_name.lower() == "dist":
            return "DISTRITO"
        elif col_name.lower() == "dist_cod":
            return "DISTRITO CODIGO"
        elif col_name.lower() == "circ":
            return "CIRCUNSCRIPCION ELECTORAL"
        elif col_name.lower() == "circ_cod":
            return "CIRCUNSCRIPCION ELECTORAL CODIGO"
        elif col_name.lower() == "local_cod":
            return "LOCAL CODIGO"
        else:
            return unidecode(col_name).upper().strip()

    # election_name = "2023_05_CCG"
    # CIRCUNSCRIPCION SENATORIAL : son 16 (se usa para eleccion de sandores)
    # DISTRITO : son 28 (se usan para eleccion de diputados)
    # CIRCUNSCRIPCION ELECTORAL : son 346 en la eleccion de convencioanles del 2023

    # llave mesa normalizada
    llave_mesa = ['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'LOCAL', 'MESA']

    # votantes
    # read excel file
    if verbose: print("\n######### Leer votantes #########")
    electors = pd.read_excel(f'{election_name}/data/{election_name}.xlsx', skiprows=6, sheet_name = 'Votantes efectivos en Chile')

    # normalize encoding of the columns content
    if verbose: print("\n### Normalizar encoding votantes ###")
    start = time.time()
    for col in electors.columns:
        if electors[col].dtype == object:
            print("\t", col)
            electors[col] = electors[col].fillna('').apply(lambda x: re.sub(' +', ' ', unidecode(x).upper().strip().replace('"','')))
    if verbose: print(f"Tiempo normalización: {time.time() - start} segundos")
    #.decode('utf-8')

    # normalize encoding of the columns names
    electors = electors.rename(columns={col: unidecode(col).upper().strip() for col in electors.columns})

    # get groups

    if verbose: print("\n### Obtener grupos ###")

    start = time.time()
    electors = electors[~((electors['RANGO ETARIO'] == '') & (electors['VOTANTES'] == 0))]
    
    # WE ARE CONSIDERING ONLY AGE GROUPS NOW
    electors['GRUPO'] = electors['RANGO ETARIO'] # + ' ' + electors['SEXO'].str[0] # full name of groups
    grupos = list(np.sort(electors['GRUPO'].unique())) # order groups: first sex then age
    electors = electors.groupby(llave_mesa + ['GRUPO']).sum().reset_index() # group by GROUPS, there may be repeated rows (NACIONALIDAD)
    electors = electors.pivot(index = llave_mesa, columns = 'GRUPO', values='VOTANTES').reset_index() # long to wide
    electors = electors.reindex(llave_mesa + grupos, axis='columns') # reindex with key 
    electors[grupos] =  electors[grupos].fillna(0).astype(int) # fill na for groups who didnt vote in that ballot box
    
    # save file with electors
    electors.to_csv(f'{election_name}/output/{election_name}_ELECTORES.csv', index=False)
    # guardar grupos
    # lista_grupos = list(natsorted(grupos))
    # with open(f'{election_name}/output/GRUPOS.pickle', "wb") as handle:
    #     pickle.dump(lista_grupos, handle)

    if verbose: print('\n### Ejemplo electores ###')
    if verbose: print(electors.head(5).to_string(index=False))
    if verbose: print(f'Tiempo obtener grupos: {time.time() - start} segundos')
    return electors

# función que dada la matriz de votantes (columnas con grupo demográficos y filas con mesas )
# simula una asignación aleatoria de votantes a mesas
def simulated_xi_sq_statistic(votantes):
    n_votantes_grupo = votantes.sum(axis = 0)
    n_votantes_mesas = votantes.sum(axis = 1)
    n_votantes_mesas_acum = np.cumsum(n_votantes_mesas)
    fila_votantes = np.concatenate([np.repeat(i, n_votantes_grupo[i]) for i in range(len(n_votantes_grupo))])
    # shuffle elements of the array
    np.random.shuffle(fila_votantes)
    fila_votantes
    electors_district_simulated = pd.DataFrame(np.zeros(votantes.shape, dtype=int), columns = votantes.columns)
    electors_district_simulated.shape
    for i in range(votantes.shape[0]):
        indice_desde = n_votantes_mesas_acum[i-1] if i > 0 else 0
        indice_hasta = n_votantes_mesas_acum[i]
        electors_district_simulated.iloc[i, :] = np.bincount(fila_votantes[indice_desde:indice_hasta], minlength = votantes.shape[1])
    chi2_stat, p_val, dof, expected = chi2_contingency(electors_district_simulated)
    return p_val

# hace test chi cuadrado con matriz de contingencia para cada distrito
def test_xi_square_per_district():
    # electors = count_groups_per_ballot("2021_11_Presidencial")
    electors = pd.read_csv("2021_11_Presidencial/output/2021_11_Presidencial_ELECTORES.csv")
    electors.head(5)
    electors.columns
    districts = electors['CIRCUNSCRIPCION ELECTORAL'].unique()
    # p_val_dict = {}
    df_p_val = pd.DataFrame(columns = ['REGION', 'DISTRITO', 'NUM_MESAS', 'P_VAL'])
    for district in districts:
        print(district)
        region = electors[electors['CIRCUNSCRIPCION ELECTORAL'] == district]['REGION'].iloc[0]
        electors_district = electors[electors['CIRCUNSCRIPCION ELECTORAL'] == district]
        electors_district = electors_district.drop(['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'LOCAL'], axis = 1)
        electors_district = electors_district.set_index('MESA')
        #expected_per_group = electors_district.sum(axis = 0) / len(electors_district)
        #if min(electors_district.sum(axis = 0) / len(electors_district))
        # Perform chi-square test
        if min(electors_district.sum(axis = 0)) == 0:
            print(electors_district)
            electors_district = electors_district.drop(electors_district.columns[electors_district.sum(axis = 0) == 0], axis = 1)
            print(electors_district)
            # print("Warning: Expected frequency less than 5"
        chi2_stat, p_val, dof, expected = chi2_contingency(electors_district)
        #p_val_dict[district] = p_val
        #df_p_val = df_p_val.append({'REGION': region, 'DISTRITO': district, 'P_VAL': p_val}, ignore_index=True)
        df_p_val = pd.concat([df_p_val, pd.DataFrame({'REGION': region, 'DISTRITO': district, 'NUM_MESAS' : len(electors_district), 'P_VAL': p_val}, index=[0])], ignore_index=True)
        # do a sampling of the multinomial
        print("P-value:", p_val)
        # p_val = simulated_xi_sq_statistic(electors_district)
        # print("P-value simulated:", p_val)
    
    df_p_val.to_csv("2021_11_Presidencial/output/2021_11_Presidencial_GROUPS_HET_P_VAL.csv", index=False)
    return df_p_val

# plot for all electoral district p-values
def plot_chi_square_pvalues(df_p_val):
    x = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    y = [100* df_p_val['P_VAL'].apply(lambda x: x < cut).mean() for cut in x]
    y2 = [100* df_p_val['P_VAL'].apply(lambda xx: ((xx >= cut / 10) or (cut == x[0])) and (xx < cut)).mean() for cut in x]
    
    # create dataframe
    df = pd.DataFrame({'x': x, 'y': y})
    print(df)

    # Create the plot
    plt.figure(figsize=(8, 4))
    plt.plot(x, y, marker='o', label='Line 1')
    plt.ylim(80, 100)  # Set the y-axis limits from 0 to 100
    plt.xscale('log')
    plt.xlabel('p-values')
    plt.ylabel('Cumulative Histogram of Districts [%]')
    plt.gca().set_yticklabels(['{:.0f}%'.format(x) for x in plt.gca().get_yticks()])  # Format y-axis as percentage
    # plt.xticks(x, [str(i) for i in x])
    plt.xticks(x)
    # save plot on hardrive as pdf
    plt.savefig("images/2021_11_Presidencial_GROUPS_HET_P_VAL.pdf")
    plt.show()

    plt.figure()
    plt.plot(x, y2, marker='o', label='Line 1')
    plt.ylim(0, 100)  # Set the y-axis limits from 0 to 100
    plt.xscale('log')
    plt.xlabel('p-values')
    plt.ylabel('Histogram of districts p-value [%]')
    plt.gca().set_yticklabels(['{:.0f}%'.format(x) for x in plt.gca().get_yticks()])  # Format y-axis as percentage
    # plt.xticks(x, [str(i) for i in x])
    plt.xticks(x)
    plt.savefig("images/2021_11_Presidencial_GROUPS_HET_P_VAL_2.pdf")
    plt.show()
    
    
    # Add grid and legend
    #plt.grid(True)
    #plt.legend()

    # Set x-axis to logarithmic scale
    #plt.xscale('log')
    
    
    # Add labels and title
    
    # Show plot
    #plt.show()


if __name__ == '__main__':
    df_p_val = test_xi_square_per_district()
    plot_chi_square_pvalues(df_p_val)


# /opt/homebrew/bin/python3.9 group_heterogeinity.py
