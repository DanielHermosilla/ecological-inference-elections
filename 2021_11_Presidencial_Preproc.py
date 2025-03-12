# pre-process election data
import pickle
import pandas as pd
import time
import numpy as np
import os   
from unidecode import unidecode
from natsort import natsorted 
import re

import sys
from p_val_mult import compute_pvalue_pickle


from EM_mult import EM_mult
from group_opt import new_group_matrix, optimze_age_groups_v2, group_election_routine
from helper_functions import round_sum_n

verbose = True # print outputs

# function that preprocess the data
def pre_process_EM(election_name, groups = ""):
    # normalized key
    llave_mesa = ['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'LOCAL', 'MESA']
    if verbose: print("\n######### Leer votos #########")
    votes = pd.read_excel(f'{election_name}/data/{election_name}.xlsx', sheet_name='Votación en Chile', skiprows = 6)
    
    # normalize encoding of the columns content
    if verbose: print("\n######### Normalizar encoding votos #########")
    start = time.time()
    for col in votes.columns:
        if votes[col].dtype == object:
            if verbose: print("\t", col)
            # Format names for consistency
            votes[col] = votes[col].fillna('').apply(lambda x: re.sub(' +', ' ', unidecode(x).upper().strip().replace('"', '').replace('.', '')))

    if verbose: print(f"Tiempo normalización: {time.time() - start} segundos")
    votes = votes.rename(columns={col: unidecode(col).upper().strip() for col in votes.columns})
    
    # create a column with the full name of the candidate
    votes['CANDIDATO'] = votes['NOMBRES'] + ' ' + votes['PRIMER APELLIDO']
    # only add CANDIDATO if PRIMER APELLIDO is not empty
    votes.loc[votes['PRIMER APELLIDO'] == '', 'CANDIDATO'] = votes['NOMBRES']
    candidatos = natsorted(votes['CANDIDATO'].unique())
    # PIVOT
    votes = votes.pivot(index=llave_mesa, columns='CANDIDATO', values='VOTOS').reset_index()
    votes = votes.rename(columns={col: unidecode(col).upper().strip() for col in votes.columns})
    votes['NULO BLANCO'] = votes['VOTOS NULOS'] + votes['VOTOS EN BLANCO'] # unir votos nulos y votos en blanco
    candidatos = candidatos + ['NULO BLANCO']
    candidatos = [c for c in candidatos if c not in ['VOTOS NULOS', 'VOTOS EN BLANCO']] # no considerar votos nulos ni en blanco, solo la suma

    # save file with votes
    # create the folder iof it does not exist
    if not os.path.exists(f'{election_name}/output{groups}'):
        os.makedirs(f'{election_name}/output{groups}')
    votes.to_csv(f'{election_name}/output{groups}/{election_name}_VOTOS.csv', index=False)
    
    if verbose: print(votes.head())

    # save candidates
    lista_candidatos = natsorted(list(candidatos))
    with open(f'{election_name}/output{groups}/CANDIDATOS.pickle', "wb") as handle:
        pickle.dump(lista_candidatos, handle)


    # voters
    # read excel file
    if verbose: print("\n######### Leer votantes #########")
    electors = pd.read_excel(f'{election_name}/data/{election_name}.xlsx', skiprows=6, sheet_name = 'Votantes efectivos en Chile')

    # normalize encoding of the columns content
    if verbose: print("\n### Normalizar encoding votantes ###")
    start = time.time()
    for col in electors.columns:
        if electors[col].dtype == object:
            print("\t", col)
            # Format names for consistency
            electors[col] = electors[col].fillna('').apply(lambda x: re.sub(' +', ' ', unidecode(x).upper().strip().replace('"', '').replace('.', '')))
    if verbose: print(f"Tiempo normalización: {time.time() - start} segundos")

    # normalize encoding of the columns names
    electors = electors.rename(columns={col: unidecode(col).upper().strip() for col in electors.columns})

    # get groups
    if verbose: print("\n### Obtener grupos ###")

    start = time.time()
    # remove rows with empty age range or 0 voters
    electors = electors[~((electors['RANGO ETARIO'] == '') & (electors['VOTANTES'] == 0))]
    
    electors['GRUPO'] = electors['RANGO ETARIO'] # + ' ' + electors['SEXO'].str[0] # full name of groups
    if groups == '_HM':
        electors['GRUPO'] = electors['GRUPO'] + ' ' + electors['SEXO'].str[0] # full name of groups
    grupos = list(np.sort(electors['GRUPO'].unique())) # order groups: first sex then age
    if groups == '_HM':
        grupos = [g for g in grupos if 'H' in g] + [g for g in grupos if 'M' in g] 
    electors = electors.groupby(llave_mesa + ['GRUPO']).sum().reset_index() # group by GROUPS, there may be repeated rows (NACIONALIDAD)
    electors = electors.pivot(index = llave_mesa, columns = 'GRUPO', values='VOTANTES').reset_index() # long to wide
    electors = electors.reindex(llave_mesa + grupos, axis='columns') # reindex with key 
    electors[grupos] =  electors[grupos].fillna(0).astype(int) # fill na for groups who didnt vote in that ballot box
    
    # IF YOU WANT TO GROUP BEFORE
    # electors_grouped, group_names_agg = new_group_matrix(electors[grupos].to_numpy(), group_combination) # get aggregated group names and aggregated matrix
    # electors[group_names_agg] = electors_grouped # assign aggregated groups


    # save file with electors
    electors.to_csv(f'{election_name}/output{groups}/{election_name}_ELECTORES.csv', index=False)
    lista_grupos = grupos.copy()
    with open(f'{election_name}/output{groups}/GRUPOS.pickle', "wb") as handle:
        pickle.dump(lista_grupos, handle)
    handle.close()

    if verbose: print('\n### Ejemplo electores ###')
    print(electors.head(5).to_string(index=False))
    if verbose: print(f'Tiempo obtener grupos: {time.time() - start} segundos')

# function that runs the EM algorithm
def run_EM(election_name, groups = ""):
    # groups : can take the values:
    #          "" : optimize groups
    #          "_40"  : group by less and greater than 40
    #          "_HM" : group by men and women
    llave_mesa = ['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'LOCAL', 'MESA']
    
    votes = pd.read_csv(f'{election_name}/output{groups}/{election_name}_VOTOS.csv')
    electors = pd.read_csv(f'{election_name}/output{groups}/{election_name}_ELECTORES.csv')
    
    # read candidatos
    with open(f'{election_name}/output{groups}/CANDIDATOS.pickle', "rb") as handle:
        candidatos = pickle.load(handle)
    handle.close()

    # read grupos
    with open(f'{election_name}/output{groups}/GRUPOS.pickle', "rb") as handle:
        grupos = pickle.load(handle)
    handle.close()

    if verbose: print("\n######### Ejecutar EM #########")

    # unidad maxima de los mismos candidatos
    nivel_agregacion_candidatos = 'REGION'
    niveles_agregacion_candidatos = votes[nivel_agregacion_candidatos].unique() # ¿ MAS DE UNA DIMENSION?

    nivel_agregacion_EM = ['CIRCUNSCRIPCION ELECTORAL']

    dict_dfs_distritos = {}

    dict_input_pvalue = {}

    if groups == "": # optimize groups
        group_combination_fixed = None
    elif groups == "_40": # group by less and greater than 40
        group_combination_fixed = [[0,1,2],[3,4,5,6,7]] # HM <= 39 y HM >= 40
    elif groups == "_HM": # group by men and women
        group_combination_fixed = [[0,1,2,3,4,5,6,7],[8,9,10,11,12,13,14,15]]
    else:
        raise ValueError("Groups not recognized")
    
    ### group_type for names
    group_type = 'age_groups'
    if groups == '_HM': 
        group_type = 'age_sex_groups'

    failed_circs = []

    # create a consolidates dataframe for the country only if groups are fixed
    if groups != "":
        df_pais_prob = pd.DataFrame()
    
    for d in niveles_agregacion_candidatos:

        # save file for distrito
        path_distrito = f'{election_name}/output{groups}/{d}'
        if not os.path.exists(path_distrito):
            os.makedirs(path_distrito)
        
        sub_votes = votes[votes[nivel_agregacion_candidatos] == d].copy()

        niveles_agregacion_EM = sub_votes[nivel_agregacion_EM].drop_duplicates().values.tolist()

        sub_electors = electors[electors['REGION'] == d].copy()
        df_distrito = sub_votes.merge(sub_electors[llave_mesa + grupos], on = llave_mesa, how = 'inner')    
        
        I = len(candidatos)

        # number of votes
        df_distrito['NUM VOTOS'] = df_distrito[candidatos].sum(axis=1)

        # ballot-boxes need at least 1 vote
        df_distrito = df_distrito[df_distrito['NUM VOTOS'] > 0]

        df_distrito['NUM MESAS'] = -1

        df_distrito['P-VALOR'] = -1
        df_distrito['LOG P-VALOR'] = -1

        for i in range(I):
            df_distrito.loc[:,f'E_{candidatos[i]}'] = -1

        for l in niveles_agregacion_EM:
            path_circ = f'{election_name}/output{groups}/{d}/{l[0]}'
            if not os.path.exists(path_circ):
                os.makedirs(path_circ)

            if verbose: print(f'\t{l[0]}')
            index_local = (df_distrito['CIRCUNSCRIPCION ELECTORAL'] == l[0])
            df_local = df_distrito.loc[index_local].copy()
            if len(df_local) == 0: # if there's no ballot-boxes
                continue

            x = df_local[candidatos].to_numpy()
            Js = df_local['NUM VOTOS'].to_numpy()
            b_ = df_local[grupos].to_numpy()

            M = x.shape[0]
            try:
                if group_combination_fixed is None:
                    if M > 1:
                        group_combinations, _ = group_election_routine(x, b_, std_thres=0.05, threshold='gc', S = 50) 
                    else: 
                        group_combinations = [[0,1,2,3,4,5,6,7]]
                else: 
                    group_combinations = group_combination_fixed.copy()
            except:
                group_combinations = group_combination_fixed.copy()
                failed_circs.append((l[0], M))
            b, group_names_agg = new_group_matrix(b_, group_combinations, group_type = group_type) # get aggregated group names and aggregated matrix

            df_local[group_names_agg] = b # assign aggregated groups

            G = len(group_names_agg) # number of groups
            
            p, iterations, t = EM_mult(x, b, verbose = False, max_iterations = 1000)
            
            for g in range(G):
                for i in range(I):
                    df_local[f'P_{candidatos[i]}_{group_names_agg[g]}'] = p[g,i]

            df_p_local = pd.DataFrame(p, columns = candidatos, index = group_names_agg)
            df_p_local.to_csv(os.path.join(path_circ,f'P_{l[0]}.csv'))
            
            r = b @ p / np.sum(b, axis = 1)[...,None] #m,i
            E_l = r * Js[...,None]

            if ~np.any(np.isnan(p)):
                for m in range(M):
                    E_l[m] = round_sum_n(E_l[m])
                df_local.loc[:,[f'E_{candidatos[i]}' for i in range(I)]] = E_l
                df_local.loc[:,[f'DIF_{candidatos[i]}' for i in range(I)]] = x - E_l
                
            df_local['NUM MESAS'] = M
            df_distrito.loc[index_local, 'NUM MESAS'] = M

            df_local.to_csv(os.path.join(path_circ,f'{l[0]}.csv'), index=False)
            if groups != "":
                df_pais_prob = pd.concat([df_pais_prob, df_local])

            # save GRUPOS 
            grupos_local = list(natsorted(group_names_agg))
            with open(os.path.join(path_circ, 'GRUPOS.pickle'), "wb") as handle:
                pickle.dump(grupos_local, handle)

            for m in range(M):

                r_m = (p.T @ b[m]) / np.sum(b[m])
                mesa = df_local.iloc[m]['MESA']
                dict_input_pvalue[l[0], mesa] = { 'r': r_m,
                                        'J': np.sum(x[m]),
                                        'x': x[m]}

        df_distrito.to_csv(os.path.join(path_distrito,f'{d}.csv'), index=False)
        dict_dfs_distritos[d] = df_distrito

    print(f'Failed circs: {failed_circs}')
    
    if groups != "":
        # save df_pais_prob to csv
        df_pais_prob.to_csv(f'{election_name}/output{groups}/{election_name}_Probs.csv', index=False)
    else:
        file_dict_dfs_distritos = f'{election_name}/output{groups}/{election_name}_dfs_distritos.pickle'
        with open(file_dict_dfs_distritos, 'wb') as handle:
            pickle.dump(dict_dfs_distritos, handle, protocol=pickle.HIGHEST_PROTOCOL)
        handle.close()

        file_dict_input_pvalues = f'{election_name}/output{groups}/{election_name}_input_pvalues.pickle'
        with open(file_dict_input_pvalues, 'wb') as handle:
            pickle.dump(dict_input_pvalue, handle, protocol=pickle.HIGHEST_PROTOCOL)
        handle.close()

    return 0

# function that computes the p-values
def run_pvalue(election_name, groups = ""):
    print("\n######### Compute p-values #########")
    file_dict_input_pvalues = f'{election_name}/output{groups}/{election_name}_input_pvalues.pickle'
    file_dict_output_pvalues = f'{election_name}/output{groups}/{election_name}_output_pvalues.pickle'
    compute_pvalue_pickle(file_in = file_dict_input_pvalues, file_out = file_dict_output_pvalues, load_bar=True, 
                          S_min = 2, S_max = 8, seed = 123)

# function that adds the p-values already computed to the dataframes
def add_pvalues(election_name, groups = ""):
    print("\n######### creating output files #########")
    
    # read output pvalues
    file_dict_output_pvalues = f'{election_name}/output{groups}/{election_name}_output_pvalues.pickle'
    with open(file_dict_output_pvalues, 'rb') as handle:
        p_value_pickle = pickle.load(handle)
    handle.close()

    # read dataframe of districts from EM
    file_dict_dfs_distritos = f'{election_name}/output{groups}/{election_name}_dfs_distritos.pickle'
    with open(file_dict_dfs_distritos, 'rb') as handle:
        dict_dfs_distritos = pickle.load(handle)
    handle.close()

    if verbose: print("#### Creando .csv para cada distrito y país. ####")
    llave_mesa = ['REGION', 'CIRCUNSCRIPCION ELECTORAL', 'LOCAL', 'MESA']
    columnas_complementarios = []
    df_pais_list = []

    lista_distritos = sorted(list(dict_dfs_distritos.keys()))
    with open(f'{election_name}/output{groups}/REGIONES.txt', "w") as output:
        output.write('\n'.join(map(str, lista_distritos)))
    output.close()

    for d in lista_distritos:
        # d = 'DISTRITO 2'
        # df_distrito = dict_dfs_distritos[d]
        df_distrito = dict_dfs_distritos[d].copy()
        for index_m, df_distrito_row in df_distrito.iterrows():
            m = df_distrito_row['MESA']
            c = df_distrito_row['CIRCUNSCRIPCION ELECTORAL']
            df_distrito.loc[index_m, 'P-VALOR'] = p_value_pickle[(c,m)]['p_value']
        df_distrito['LOG P-VALOR'] = np.log10(df_distrito['P-VALOR'])
        df_distrito.to_csv(f'{election_name}/output{groups}/{d}/{d}.csv', index=False)
        df_pais_list.append(df_distrito[llave_mesa + ['NUM MESAS', 'NUM VOTOS'] + columnas_complementarios + ['P-VALOR', 'LOG P-VALOR']])
        
        # csv por local
        circ_list = list(df_distrito['CIRCUNSCRIPCION ELECTORAL'].unique())

        # guardar circunscripciones
        with open(f'{election_name}/output{groups}/{d}/CIRCUNSCRIPCIONES.txt', "w") as output:
            output.write('\n'.join(map(str, circ_list)))
        output.close()

    df_pais = pd.concat(df_pais_list)
    df_pais.to_csv(f'{election_name}/output{groups}/{election_name}_PAIS.csv', index=False)
    
    return None

# run with group optimization
def case_group_optimization():
    groups = ""
    election_name = "2021_11_Presidencial"
    start_total = time.time()
    
    if verbose: print('PROCESS EM')
    t_start = time.time()
    pre_process_EM(election_name, groups = groups)
    print('Time pre_process_EM: ', time.time() - t_start, ' seconds')


    if verbose: print('RUN EM')
    t_start = time.time()
    run_EM(election_name, groups = groups)
    print('Time run_EM: ', time.time() - t_start, ' seconds')

    if verbose: print('RUN P-VAL')
    t_start = time.time()
    run_pvalue(election_name, groups = groups)
    print('Time run_pvalue: ', time.time() - t_start, ' seconds')

    if verbose: print('ADD P-VAL')
    t_start = time.time()
    add_pvalues(election_name, groups = groups)
    print('Time add_pvalues: ', time.time() - t_start, ' seconds')

    if verbose: print('TOTAL TIME: ', time.time() - start_total, ' seconds')

# run for two groups: <40 and ≥40
def case_group_age():
    groups = "_40"
    election_name = "2021_11_Presidencial"
    
    start_total = time.time()
    
    if verbose: print('PROCESS EM')
    t_start = time.time()
    pre_process_EM(election_name, groups = groups)
    print('Time pre_process_EM: ', time.time() - t_start, ' seconds')

    if verbose: print('RUN EM')
    t_start = time.time()
    run_EM(election_name, groups = groups)
    print('Time run_EM: ', time.time() - t_start, ' seconds')

    if verbose: print('TOTAL TIME: ', time.time() - start_total, ' seconds')

# run for two groups: male and female
def case_group_sex():
    groups = "_HM"
    election_name = "2021_11_Presidencial"
    
    start_total = time.time()

    if verbose: print('PROCESS EM')
    t_start = time.time()
    pre_process_EM(election_name, groups = groups)
    print('Time pre_process_EM: ', time.time() - t_start, ' seconds')

    if verbose: print('RUN EM')
    t_start = time.time()
    run_EM(election_name, groups = groups)
    print('Time run_EM: ', time.time() - t_start, ' seconds')

    if verbose: print('TOTAL TIME: ', time.time() - start_total, ' seconds')

if __name__ == '__main__':
    case = "" if len(sys.argv) == 1 else sys.argv[1]
    if case == "":
        case_group_optimization()
    elif case == "age":
        case_group_age()
    elif case == "sex":
        case_group_sex()
    else:
        print("Case not recognized")
        exit(1)

# python3 2021_11_Presidencial_Preproc.py
# python3 2021_11_Presidencial_Preproc.py age
# python3 2021_11_Presidencial_Preproc.py sex

# 2445.9 seconds for G3
# 600 seconds with parallel

# 4719 segs