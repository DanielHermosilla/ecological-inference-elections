import pickle
import numpy as np
import pandas as pd
import sys

# function that read the instances results and return a dataframe
def read_results():

    with_replace = False

    J_list = [100] # personas
    M_list = [50] # mesas
    # G_list = [2,3,4] # grupos
    # I_list = [2,3,4,5,10] # candidatos 
    G_list = [2]
    I_list = [2]
    lambda_list = [0.5]
    # seed_list = [i+1 for i in range(20)]
    seed_list = [1]
    
    instances = []
    n_instances = len(J_list)*len(M_list)*len(G_list)*len(I_list)*len(seed_list)
    # EM_method_names = ["full", "simulate_100", "simulate_1000", "cdf", "pdf", "mult"]
    if with_replace:
        EM_method_names = ["full", "simulate_100_with_replace", "simulate_1000_with_replace", "cdf", "pdf", "mult"]
    else:
        EM_method_names = ["full", "simulate_100", "simulate_1000", "cdf", "pdf", "mult"]

    df_list = []
    for J in J_list:
        for M in M_list:
            for G in G_list:
                for I in I_list:
                    for lambda_ in lambda_list:
                            for method in EM_method_names:
                                for seed in seed_list:
                                    if 'simulate' in method:
                                        # J = 100
                                        # M = 50
                                        # G = 2
                                        # I = 2
                                        # lambda_ = 0.5
                                        # method = 'cdf'
                                        # seed = 1
                                        
                                        pickle_path = f'results/J{J}_M{M}_G{G}_I{I}_lambda{int(100*lambda_)}/{method}_cv1000/{seed}.pickle'
                                        # pickle_path = f'/Users/charlesthraves/Dropbox (MIT)/Research/SERVEL/Mesas Outlier 2/Codigos/Pablo/elections-model/results/J{J}_M{M}_G{G}_I{I}_lambda{int(100*lambda_)}/{method}_cv1000/{seed}.pickle'
                                        # check if path exists
                                        # try:
                                        with open(pickle_path, 'rb') as f:
                                            data = pickle.load(f)
                                        
                                        p = data['p']
                                
                                        p_est = data['p_est']
                                        mean_error = np.mean(np.abs(p-p_est))
                                        max_error = np.max(np.abs(p-p_est))
                                        time = data['time']
                                    
                                        if 'simulate' in method:
                                            simulate_time = data['simulation_time']
                                        else:
                                            simulate_time = np.nan
                                        iter = data['iterations']
                                        end = data['end']
                                        df_list.append([J,M,G,I,lambda_,method,seed,time,simulate_time,iter,end,mean_error,max_error])
                                    # except:
                                        
                                    #     print(f'Error in {pickle_path}')
                                    #     continue
    
    df = pd.DataFrame(df_list, columns=['J','M','G','I','lambda','method','seed','time','simulate_time','iter','end','mean_error','max_error'])                                  
    return df

# function that prints the running time table
def latex_table_time():

    # choose method order
    with_replace = False
    with_replace_str = '_with_replace' if with_replace else ''
    method_order = ['full', f'simulate_100{with_replace_str}', f'simulate_1000{with_replace_str}', 'cdf', 'pdf', 'mult']
    
    # decimal digits
    decimals_method = [1, 1, 1, 2, 2, 4]

    # get the data of results
    df = read_results()

    # replace time for nan if end = -1 
    df.loc[df['end'] == 0, 'time'] = np.nan
    df['total_time'] = df['time'] + df['simulate_time'].fillna(0)

    ### include the simulation time ##
    # grab the average time with J,M,G,I as row and method as column
    df_time = df.groupby(['J','M','G','I','method']).mean().reset_index().pivot(index=['J','M','G','I'], columns='method', values='total_time').reset_index()
    df_time[['J', 'M', 'G', 'I']] = df_time[['J', 'M', 'G', 'I']].astype(int)

    # round to 4 decimals for methods
    # df_time[method_order] = df_time[method_order].round(3)
    # round to decimals_method respectively
    df_time[method_order] = df_time[method_order].apply(lambda x: np.round(x, decimals_method[method_order.index(x.name)]))

    for method, decimals in zip(method_order, decimals_method):
        df_time[method] = df_time[method].apply(lambda x: f'{x:.{decimals}f}')

    # df_time[method_order] = df_time[method_order].applymap(lambda x: "{:.3f}".format(x))

    # add separator for thousands and put between brackets for those that have a comma
    # df_time[method_order] = df_time[method_order].applymap(lambda x: f'{float(x):,.3f}' if float(x) >= 1000 else x) #.replace(',', ' ')

    # display floats with 3 decimals even if they are 0'seed

    # put numbers with commas between brackets
    # df_time[method_order] = df_time[method_order].applymap(lambda x: f'{{{x}}}' if ',' in x else x)

    # replace 'nan' with '-' 
    df_time = df_time.replace('nan', '{-}')



    # reorder columns
    J = 100
    df_time_latex = pd.DataFrame()
    # change order of columns and use english notation
    df_time_latex[['C', 'G'] + method_order] = df_time.loc[df_time['J'] == J, ['I', 'G'] + method_order]
    # sort by C and G
    df_time_latex = df_time_latex.sort_values(['C', 'G'])


    df_time_latex.loc[df_time_latex['C'].duplicated(),['C']] = ''

    # add empty columns between 'C' and methods
    df_time_latex.insert(2, 'empty1', '')

    # transform to latex
    print(df_time_latex.to_latex(index=False, escape=False, column_format='cccccccccc', float_format=":.%4f"))

    return 0

def latex_table_simulation_time():
    method_order_sim = ['simulate_timesimulate_100_with_replace', 'timesimulate_100_with_replace', 'simulate_timesimulate_1000_with_replace', 'timesimulate_1000_with_replace']
    df = read_results()

    # grab the average time with J,M,G,I as row and method as column
    df_sim = df.loc[df['method'].str.contains('simulate_1')].groupby(['J','M','G','I','method']).mean().reset_index().pivot(index=['J','M','G','I'], columns='method', values=['time', 'simulate_time']).reset_index()

    df_sim[['J', 'M', 'G', 'I']] = df_sim[['J', 'M', 'G', 'I']].astype(int)
    df_sim.columns = [f'{col[0]}{col[1]}' for col in df_sim.columns]

    # round to 4 decimals for methods
    df_sim[method_order_sim] = df_sim[method_order_sim].round(3)
    df_sim[method_order_sim] = df_sim[method_order_sim].applymap(lambda x: "{:.3f}".format(x))

    # add separator for thousands and put between brackets for those that have a comma
    df_sim[method_order_sim] = df_sim[method_order_sim].applymap(lambda x: f'{float(x):,.3f}' if float(x) >= 1000 else x) #.replace(',', ' ')

    # display floats with 3 decimals even if they are 0'seed

    # put numbers with commas between brackets
    # df_sim[method_order] = df_sim[method_order].applymap(lambda x: f'{{{x}}}' if ',' in x else x)

    # replace 'nan' with '-' 
    df_sim = df_sim.replace('nan', '{-}')

    # reorder columns
    J = 100
    df_time_sim_latex = pd.DataFrame()

    # change order of columns and use english notation
    df_time_sim_latex[['C', 'G'] + method_order_sim] = df_sim.loc[df_sim['J'] == J, ['I', 'G'] + method_order_sim]
    # sort by C and G
    df_time_sim_latex = df_time_sim_latex.sort_values(['C', 'G'])


    df_time_sim_latex.loc[df_time_sim_latex['C'].duplicated(),['C']] = ''

    # add empty columns between 'C' and methods
    df_time_sim_latex.insert(2, 'empty1', '')

    # transform to latex
    print(df_time_sim_latex.to_latex(index=False, escape=False, float_format=":.%4f"))
    return 0

def latex_table_error():
    with_replace = False
    with_replace_str = '_with_replace' if with_replace else ''
    method_order = ['full', f'simulate_100{with_replace_str}', f'simulate_1000{with_replace_str}', 'cdf', 'pdf', 'mult']
    
    df = read_results()
    
    # check this case for the mean error (when it did'nt finish)
    df.loc[df['end'] == 0,'mean_error'] = np.nan

    # CG   

    # grab the average time with J,M,G,I as row and method as column
    df_error = df.groupby(['J','M','G','I','method']).mean().reset_index().pivot(index=['J','M','G','I'], columns='method', values='mean_error').reset_index()
    df_error[['J', 'M', 'G', 'I']] = df_error[['J', 'M', 'G', 'I']].astype(int)
    # 
    # round to 4 decimals for methods
    df_error[method_order] = df_error[method_order].round(3)
    df_error[method_order] = df_error[method_order].applymap(lambda x: "{:.3f}".format(x))

    # add separator for thousands and put between brackets for those that have a comma
    df_error[method_order] = df_error[method_order].applymap(lambda x: f'{float(x):,.3f}' if float(x) >= 1000 else x) #.replace(',', ' ')

    # display floats with 3 decimals even if they are 0'seed

    # put numbers with commas between brackets
    # df_error[method_order] = df_error[method_order].applymap(lambda x: f'{{{x}}}' if ',' in x else x)

    # replace 'nan' with '-' 
    df_error = df_error.replace('nan', '{-}')



    # reorder columns
    J = 100
    df_error_latex = pd.DataFrame()
    # change order of columns and use english notation
    df_error_latex[['C', 'G'] + method_order] = df_error.loc[df_error['J'] == J, ['I', 'G'] + method_order]
    # sort by C and G
    df_error_latex = df_error_latex.sort_values(['C', 'G'])


    df_error_latex.loc[df_error_latex['C'].duplicated(),['C']] = ''

    # add empty columns between 'C' and methods
    df_error_latex.insert(2, 'empty1', '')

    # transform to latex
    print(df_error_latex.to_latex(index=False, escape=False, column_format='cccccccccc', float_format=":.%4f"))

    # # replace I by C
    # df_error['C'] = df_time['I']

def latex_table_error_v2():
    with_replace = False
    with_replace_str = '_with_replace' if with_replace else ''
    method_order = ['full', f'simulate_100{with_replace_str}', f'simulate_1000{with_replace_str}', 'cdf', 'pdf', 'mult']
    
    df = read_results()
    # CG   

    # grab the average time with J,M,G,I as row and method as column
    df_error = df.groupby(['J','M','G','I','method']).mean().reset_index().pivot(index=['J','M','G','I'], columns='method', values='mean_error').reset_index()
    df_error[['J', 'M', 'G', 'I']] = df_error[['J', 'M', 'G', 'I']].astype(int)
    # 
    # round to 4 decimals for methods
    df_error[method_order] = df_error[method_order].round(3)
    df_error[method_order] = df_error[method_order].applymap(lambda x: "{:.3f}".format(x))

    # add separator for thousands and put between brackets for those that have a comma
    df_error[method_order] = df_error[method_order].applymap(lambda x: f'{float(x):,.3f}' if float(x) >= 1000 else x) #.replace(',', ' ')

    # display floats with 3 decimals even if they are 0'seed

    # put numbers with commas between brackets
    # df_error[method_order] = df_error[method_order].applymap(lambda x: f'{{{x}}}' if ',' in x else x)

    # replace 'nan' with '-' 
    df_error = df_error.replace('nan', '{-}')



    # reorder columns
    J = 100
    df_error_latex = pd.DataFrame()
    # change order of columns and use english notation
    df_error_latex[['G', 'C'] + method_order] = df_error.loc[df_error['J'] == 100, ['G', 'I'] + method_order]
    # sort by C and G
    df_error_latex = df_error_latex.sort_values(['G', 'C'])


    df_error_latex.loc[df_error_latex['G'].duplicated(),['G']] = ''

    # add empty columns between 'C' and methods
    df_error_latex.insert(2, 'empty1', '')

    # transform to latex
    print(df_error_latex.to_latex(index=False, escape=False, column_format='cccccccccc', float_format=":.%4f"))

    # # replace I by C
    # df_error['C'] = df_time['I']
    return 0

# main
if __name__ == '__main__':
    table_type = sys.argv[1]
    if table_type == 'time':
        latex_table_time()
    elif table_type == 'simulation_time':
        latex_table_simulation_time()
    elif table_type == 'error':
        latex_table_error()
    
# python3 analize_EM.py time
# python3 analize_EM.py simulation_time
# python3 analize_EM.py error