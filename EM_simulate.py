import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
import pickle
import os
import time
from tqdm import tqdm
from EM_algorithm import EM_algorithm, get_p_est
from instance_maker import gen_instance_v3
import json

""" Creation of samples for the Hit and Run method """


def polytope_elections(n, b):
    """ """
    I_size, G_size = len(n), len(b)

    # Ax <= c representes the inequalities of the polytope
    A = np.zeros(
        (
            G_size + I_size - 1 + (G_size - 1) * (I_size - 1),
            (G_size - 1) * (I_size - 1),
        ),
        dtype=int,
    )
    c = np.zeros(G_size + I_size - 1 + (G_size - 1) * (I_size - 1), dtype=int)

    for g in range(G_size - 1):
        # first G-1 inequalities
        ineq_index = [g * (I_size - 1) + i for i in range(I_size - 1)]
        A[g][ineq_index] = 1.0
        c[g] = b[g]
    for i in range(I_size - 1):
        # next I-1 inequalities
        ineq_index = [g * (I_size - 1) + i for g in range(G_size - 1)]
        A[G_size - 1 + i][ineq_index] = 1.0
        c[G_size - 1 + i] = n[i]
    # last inequality
    A[G_size + I_size - 2][:] = -1.0
    c[G_size + I_size - 2] = -(np.sum(n[:-1]) - b[G_size - 1])  # revisar

    for g in range(G_size - 1):
        for i in range(I_size - 1):
            A[G_size + I_size - 2 + 1 + g * (I_size - 1) + i][
                g * (I_size - 1) + i
            ] = -1.0

    return A, c


def random_direction(var_size, l, u):
    permut = np.random.permutation(var_size)
    v = np.random.randint(l, u)
    if (v == 0).all():
        return random_direction(var_size, l, u)
    return permut, v


def generate_z_m(n_m, p, b_m, G_size, I_size):
    z_m = np.zeros(shape=(G_size, I_size), dtype=int)
    b_sim_m_n = b_m.copy()  # Usar punteros al implementar en C
    votes_left_n = n_m.copy()  # Usar punteros al implementar en C
    total_votes = np.sum(n_m)
    # let's fill Z
    for i in range(total_votes):
        p_i = votes_left_n / np.sum(votes_left_n)
        choose_i = np.random.choice(I_size, p=p_i)
        p_g_cond_i = p[:, choose_i] * b_sim_m_n / np.sum(p[:, choose_i] * b_sim_m_n)
        g_vote = np.random.choice(G_size, p=p_g_cond_i)
        z_m[g_vote, choose_i] += 1
        b_sim_m_n[g_vote] -= 1
        votes_left_n[choose_i] -= 1
    return z_m


def generate_z_m_v2(n_m, b_m, G_size, I_size):
    """
    Generate the matrix of all possible outcomes.

    Parameters:
        n_m (numpy.ndarray): Total amount of votes from each ballot box.
        b_m (numpy.ndarray): Total amount of votes per demographic group for each ballot box.
        G_size (int): The amount of demographic groups.
        I_size (int): The amount of candidates.

    Returns:
        numpy.ndarray: Matrix with all the posible outcomes.
    """
    z_m = np.zeros(shape=(G_size, I_size), dtype=int)
    b_m_copy = b_m.copy()  # Usar punteros al implementar en C
    n_m_copy = n_m.copy()  # Usar punteros al implementar en C

    # generate a fast initial permutation
    for g in range(G_size):
        for i in range(I_size):
            z_m[g, i] = min(b_m_copy[g], n_m_copy[i])
            b_m_copy[g] -= z_m[g, i]
            n_m_copy[i] -= z_m[g, i]
    return z_m


# transform z_m to a reduced vector of variables
def z_to_polytope_var(z_m):
    """
    Transform the sampling vector to an array.

    Parameters:
        z_m (numpy.ndarray): The two dimensional matrix of possible results.

    Returns:
        numpy.ndarray: Flattened array


    """
    return z_m[:-1, :-1].flatten()


# transform a reduced vector of variables to z_m
def polytope_var_to_z(x_m, n_m, b_m, G_size, I_size):
    """
    Given a sample (polytope), it computes the matrix of all possible outcomes.

    Parameters:
        x_m (numpy.ndarray): Matrix of votes towards a candidate for each ballot box.
        n_m (numpy.ndarray): Total amount of votes from each ballot box. (Might be redundant with x_m).
        b_m (numpy.ndarray): Total amount of votes per demographic group for each ballot box.
        G_size (int): The amount of demographic groups.
        I_size (int): The amount of candidates.

    Returns:
        numpy.ndarray: Matrix with all the posible outcomes.

    """
    Z_m = np.zeros((G_size, I_size), dtype=int)
    Z_m[:-1, :-1] = x_m.reshape((G_size - 1, I_size - 1))
    for g in range(G_size - 1):
        Z_m[g, I_size - 1] = b_m[g] - np.sum(Z_m[g, :])
    for i in range(I_size):
        Z_m[G_size - 1, i] = n_m[i] - np.sum(Z_m[:, i])
    return Z_m


def compute_beta_z(z, log_p, lgac_n, G_size, I_size):
    """
    Compute a starting point within the matrix of possible outcomes.

    Parameters:
        z (numpy.ndarray): Matrix with all the posible outcomes.
        log_p (numpy.ndarray): (Log of probabilities of Z ¿)
        lgac_n (numpy.ndarray): (Array with normalization values ¿)
        G_size (int): The amount of demographic groups.
        I_size (int): The amount of candidates.

    Returns:
        numpy.ndarray: Array with a convenient starting point for the initialization of the H&R algorithm.
    """
    # beta_z = 0.0
    # for g in range(G_size):
    #    beta_z += np.sum([hlg_p_gin[g][i][z[g,i]] - lgac_n[z[g,i]] for i in range(I_size)])
    beta_z = np.sum(log_p * z) + np.sum(
        [-lgac_n[z[g, i]] for g in range(G_size) for i in range(I_size)]
    )
    return beta_z


# @profile
def hit_and_run_matrix(
    n_m,
    p,
    b_m,
    G_size,
    I_size,
    samples=1000,
    step_size=100,
    p_distribution=False,
    load_bar=False,
    verbose=False,
    unique=True,
    seed=123,
):
    """
    Matrix with the outcomes of the Hit and Run algorithm.

    Parameters:
        n_m (numpy.ndarray): Total amount of votes per ballot box.
        b_m (numpy.ndarray): Matrix of votes towards a candidate for each ballot box.
        G_size (int): The amount of demographic groups.
        I_size (int): The amount of candidates.
        samples (int, optional): Amount of samples (equivalent of "S" in the pseudo-code).
            default value: 1000
        step_size (int, optional): The size of the step (equivalent of "M" in the pseudo-code).
            default value: 100
        load_bar (bool, optional): Print a loading bar showing the progress.
            default value: False
        unique (bool, optional): Consider only unique appearances.
            default value: True
        seed (int, optional): The randomizing seed.
            default value: 123

        Returns:
            Returns a list of outcomes.


    """
    # z_0 = generate_z_m(n_m, p, b_m, G_size, I_size)
    # print('getting_first_point...')
    if seed is not None:
        np.random.seed(seed)

    Z_list = []
    Z = generate_z_m_v2(n_m, b_m, G_size, I_size)
    dir_samples = samples * step_size
    f1 = np.random.randint(0, G_size, size=dir_samples)  # fila casilla 1
    f2 = np.random.randint(0, G_size - 1, size=dir_samples)  # fila casilla 2
    f2[f2 >= f1] += 1
    c1 = np.random.randint(0, I_size, size=dir_samples)  # columna casilla 1
    c2 = np.random.randint(0, I_size - 1, size=dir_samples)  # columna casilla 2
    c2[c2 >= c1] += 1
    for s in tqdm(range(dir_samples), disable=not load_bar):
        z_neg_1 = Z[f1[s], c1[s]] - 1
        z_neg_2 = Z[f2[s], c2[s]] - 1
        if (z_neg_1 >= 0) and (z_neg_2 >= 0):
            Z[f1[s], c1[s]] = z_neg_1
            Z[f2[s], c2[s]] = z_neg_2
            Z[f2[s], c1[s]] += 1
            Z[f1[s], c2[s]] += 1
        # else:
        #     Z[f1[s],c1[s]] += 0
        #     Z[f2[s],c2[s]] += 0
        #     Z[f2[s],c1[s]] += 0
        #     Z[f1[s],c2[s]] += 0
        if s % step_size == 0:
            Z_list.append(Z.copy())
    # return  Z_list
    if unique:
        return np.unique(Z_list, axis=0)
    return Z_list


def simulate_Z(
    n,
    p,
    b,
    samples=1000,
    step_size=100,
    p_distribution=False,
    save=False,
    name=None,
    parallel=True,
    load_bar=False,
    verbose=False,
    unique=True,
    seed=123,
):
    np.random.seed(seed)  # set seed
    M_size, G_size, I_size = b.shape[0], b.shape[1], n.shape[1]
    # Parallelize the process
    if parallel:
        print("parallel")
        # n_m, p, b_m, G_size, I_size, samples = 100, step_size = 20, lambda_factor = 1.0, p_distribution = False, load_bar = False, verbose = False
        z_args = [
            (
                n[m],
                p,
                b[m],
                G_size,
                I_size,
                samples,
                step_size,
                p_distribution,
                load_bar,
                verbose,
                unique,
                seed * 100 + m,
            )
            for m in range(M_size)
        ]
        with Pool() as pool:
            Z = pool.starmap(hit_and_run_matrix, z_args)

    # Run sequentially
    else:
        Z = []
        for m in range(M_size):
            if verbose:
                print(f"M = {m}")
            # n_m, p, b_m, G_size, I_size, samples = 100, step_size = 20, lambda_factor = 1.0, p_distribution = False, load_bar = False, verbose = False
            Z_m = hit_and_run_matrix(
                n[m],
                p,
                b[m],
                G_size,
                I_size,
                samples,
                step_size,
                p_distribution,
                load_bar,
                verbose,
                unique,
                seed,
            )
            Z.append(Z_m)

    # Save the results
    if save:
        # Save the results in a pickle file
        if name is None:
            name = f"Z_{M_size}_{G_size}_{I_size}_{samples}"

        file_path = os.path.join(os.getcwd(), "Z_instances", f"{name}.pickle")

        with open(f"{file_path}", "wb") as f:
            pickle.dump(Z, f)

    return Z


def compute_q_simulate(
    n, p, b, Z, alfa_Z, parallel=False, parallel_2=False, parallel_3=False
):

    M_size, G_size, I_size = b.shape[0], b.shape[1], n.shape[1]
    log_p = np.log(p)

    if parallel:
        q_args = [(b[m], Z[m], I_size, G_size, alfa_Z[m], log_p) for m in range(M_size)]

        with Pool() as pool:
            q = pool.starmap(compute_qm_simulate, q_args)
        return np.array(q)

    q = np.zeros(shape=(M_size, G_size, I_size))
    for m in range(M_size):
        q[m] = compute_qm_simulate(b[m], Z[m], I_size, G_size, alfa_Z[m], log_p)

    return q


# CAMBIAR ALFA_Z
def compute_qm_simulate(b_m, Z_m, I_size, G_size, alfa_Z_m, log_p):
    probs_zm = np.array(
        [
            compute_qzm(Z_m[s], I_size, G_size, alfa_Z_m[s], log_p)
            for s in range(len(Z_m))
        ]
    )
    Z_m_pond = np.sum(Z_m * probs_zm[..., None, None], axis=0)
    cumulative_prob = np.sum(probs_zm)
    # for z_m in Z_m:
    #    prob_zm = compute_qzm(z_m, p, b_m, I_size, G_size, lgac_n, hlg_p_gin)
    #    q_m += prob_zm * z_m
    #    cumulative_prob += prob_zm

    q_m = Z_m_pond / b_m[..., None] / cumulative_prob

    q_m[np.isnan(q_m)] = 0  # REVISAR
    # if np.any(cumulative_prob) == np.nan:
    #    print("cumulative_prob = 0")
    return q_m


def compute_qzm(z, I_size, G_size, alfa_z_m, log_p):
    prob = np.sum(z * log_p) + alfa_z_m
    return np.exp(prob)


def compute_alfa_Z(Z, b, G_size, I_size):
    lgac_n = [
        sum([np.log(max(k, 1)) for k in range(j + 1)]) for j in range(np.max(b) + 1)
    ]
    alfa_Z = [
        [np.sum([lgac_n[b[m, g]] for g in range(G_size)]) for _ in Z[m]]
        for m in range(len(Z))
    ]
    for m in range(len(Z)):
        for index, z_m in enumerate(Z[m]):
            alfa_Z[m][index] -= np.sum(
                [lgac_n[z_m[g, i]] for i in range(I_size) for g in range(G_size)]
            )
            # alfa_Z[m] += np.sum([lgac_n[b[m,g]] + np.sum([lgac_n[z_m[i]] - lgac_n[z_m[i]-1] for i in range(I_size)]) for g in range(G_size)])
    return alfa_Z


def compute_p(q, b):
    num = np.sum(np.multiply(q, b[..., None]), axis=0)
    dem = np.sum(b, axis=0)[..., None]
    return num / dem


# def EM_simulate(X, b, convergence_value = 0.0001, max_iterations = 100,
#                  p_method = 'proportional', simulate = False, instance = None, load_bar = True, verbose = True):
#     M_size, I_size = X.shape
#     G_size = b.shape[1]
#     J_mean = np.round(np.mean(np.sum(X, axis = 1)),1)
#     #print('-'*100)
#     if verbose: print('M =',M_size,' G =',G_size,' I =',I_size,' J =',J_mean,' delta =',convergence_value)
#     #print('-'*100)
#     if p_method == 'uniform':
#         p_est = np.full((G_size, I_size), 1/I_size)
#     if p_method == 'proportional':
#         p_est = np.array([np.sum(X, axis = 0) / np.sum(X) for g in range(G_size)])
#         # p_est = np.array([[sum([X[m,i] for m in range(M_size)]) /
#         #                  sum([X[m] for m in range(M_size)]) for i in range(I_size)] for g in range(G_size)])
#     if p_method == 'group_proportional':
#         p_mesas = X/np.sum(X,axis=1)[...,None]
#         p_est = np.zeros((G_size, I_size))
#         for g in range(G_size):
#             for i in range(I_size):
#                 p_est[g,i] = np.sum(p_mesas[:,i]*b[:,g])/np.sum(b[:,g])

#     if simulate:
#         #print('Simulating...')
#         Z = simulate_Z(X, p_est, b, verbose = verbose, load_bar = load_bar)

#     else:
#         #start = time.time()
#         if instance == None:
#             instance = input('Cargue instancia:')
#         file_path = os.path.join(os.getcwd(), 'Z_instances', f'{instance}.pickle')
#         with open(f'{file_path}.pickle', 'rb') as f:
#             Z = pickle.load(f)
#     alfa_Z = compute_alfa_Z(Z, b, G_size, I_size)

#     if verbose:
#         print('-'*100)
#         print('EM-algorithm')
#         print('-'*100)
#     start = time.time()
#     for i in tqdm(range(1,max_iterations+1), disable = not load_bar):
#         q = compute_q_simulate(X, p_est, b, Z, alfa_Z)
#         p_new = compute_p(q,b)
#         if verbose:
#             print('-'*50)
#             print(p_est)
#             print('-'*50)
#             print('Δ: ', np.max(np.abs(p_new-p_est)))
#         #print(f'{end-start} seconds')
#         if (np.abs(p_new-p_est) < convergence_value).all():
#             end = time.time()
#             if verbose: print(f'Convergence took {i} iterations and {end-start} seconds.')
#             return p_new, i, end-start
#         p_est = p_new.copy()
#     end = time.time()
#     if verbose: print(f'Did not reach convergence after {i} iterations and {end-start} seconds. \n')
#     return p_est, i, end-start


def EM_simulate(
    X,
    b,
    p_est=None,
    convergence_value=0.0001,
    max_iterations=100,
    p_method="group_proportional",
    simulate=False,
    Z=None,
    unique=True,
    instance=None,
    samples=1000,
    step_size=100,
    load_bar=True,
    verbose=True,
    dict_results={},
    save_dict=False,
    dict_file=None,
):
    G_size, I_size = b.shape[1], X.shape[1]
    if p_est is None:
        p_est = get_p_est(X, b, p_method)

    dict_results["samples"] = samples
    dict_results["step_size"] = step_size

    if simulate:
        if verbose:
            print("Simulating...")
        start = time.time()
        Z = simulate_Z(
            X,
            p_est,
            b,
            samples=samples,
            step_size=step_size,
            verbose=verbose,
            load_bar=load_bar,
            unique=unique,
        )
        simulate_time = time.time() - start
        dict_results["simulation_time"] = simulate_time
        # dict_results['Z'] = Z

    else:
        # start = time.time()
        if Z == None:
            if instance == None:
                instance = input("Cargue instancia:")
            file_path = os.path.join(os.getcwd(), "Z_instances", f"{instance}.pickle")
            with open(f"{file_path}.pickle", "rb") as f:
                Z = pickle.load(f)
    alfa_Z = compute_alfa_Z(Z, b, G_size, I_size)

    q_method = lambda X, p, b: compute_q_simulate(X, p, b, Z, alfa_Z)

    return EM_algorithm(
        X,
        b,
        p_est,
        q_method,
        convergence_value,
        max_iterations,
        load_bar,
        verbose,
        dict_results,
        save_dict,
        dict_file,
    )


if __name__ == "__main__":
    # generate instance to test Z
    G_list = [2, 4]
    I_list = [2, 10]
    M = 50
    instances = 20
    J = 100
    step_size = 100
    samples = 10000
    for G in G_list:
        for I in I_list:
            for seed in range(instances):
                gen_instance_v3(
                    G,
                    I,
                    M,
                    J,
                    name=f"instance_G{G}_I{I}_M{M}_J{J}_seed{seed}",
                    terminar=False,
                    seed=seed,
                )

                # load instance with json instead of pickle
                with open(
                    f"instances/instance_G{G}_I{I}_M{M}_J{J}_seed{seed}.json", "r"
                ) as f:
                    data = json.load(f)
                X = np.array(data["n"])
                b = np.array(data["b"])
                p = np.array(data["p"])
                print(G, I)
                # simulate Z
                time_start = time.time()
                # Z = hit_and_run_matrix(X[0], p, b[0], G, I, samples = 1000, step_size=300, load_bar=True, verbose=True)
                Z = simulate_Z(
                    X,
                    p,
                    b,
                    samples=samples,
                    step_size=step_size,
                    verbose=True,
                    load_bar=True,
                    save=True,
                    name=f"Z_instance_G{G}_I{I}_M{M}_J{J}_step{step_size}_S{samples}_seed{seed}_sequence",
                    unique=False,
                    parallel=True,
                )
                print(f"Time elapsed: {time.time() - time_start} seconds")

        # # save as pickle
    # with open(f"Z_instance_G{G}_I{I}_M{M}_J{J}_seed{0}.pickle", 'wb') as f:
    #     pickle.dump(Z, f)

    # # Z = [Z[s].to_list() for s in range(len(Z))]
    # # # save Z as json instead of pickle
    # # with open(f"Z_instance_G{G}_I{I}_M{M}_J{J}_seed{0}.json", 'w') as f:
    # #     json.dump(Z, f)

    # load instance with json instead of pickle

    # from instance_maker import gen_instance_v3

    # G,I,M,J = 10,10,5,100

    # gen_instance_v3(G,I,M,J, name=f"instance_G{G}_I{I}_M{M}_J{J}_seed{0}", terminar=False, seed = 0)

    # with open(f"instances/instance_G{G}_I{I}_M{M}_J{J}_seed{0}.json", 'r') as f:
    #     data = json.load(f)
    # X = np.array(data["n"])
    # b = np.array(data["b"])
    # p = np.array(data["p"])
    # import time
    # time_start = time.time()
    # Z = simulate_Z(X, p, b, samples = 100, step_size = 3000, verbose=True, load_bar=True, save=False, unique=True)
    # print(f"Time elapsed: {time.time() - time_start} seconds")

    pass

