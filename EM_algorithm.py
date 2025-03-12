import numpy as np
import time
from time import perf_counter
from tqdm import tqdm
import pickle
from helper_functions import *

epsilon = 1e-10  # Small stabilizing value

"""Main EM algorithm functions to estimate probability"""


def compute_p(q, b):
    """
    Computes the optimal solution of the M-step.

    Parameters:
        q (numpy.ndarray): Matrix of dimension (bxgxc) that represents the probability that a voter of group "g" in ballot box "b" voted for candidate "c" conditional on the observed result.
        b (numpy.ndarray): Matrix of dimension (bxg) that stores the amount of votes from demographic group "g".

    Returns:
        float: The optimal probability from the M-step.
    """
    num = np.sum(np.multiply(q, b[..., None]), axis=0)
    dem = np.sum(b, axis=0)[..., None]
    return num / dem


def get_p_est(X, b, p_method):
    """
    Computes the initial probability of the EM algorithm

    Parameters:
        X (numpy.ndarray): Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
        b (numpy.ndarray): Matrix of dimension (bxg) that stores the amount of votes from demographic group "g".
        p_method (string): The method for calculating the initial parameter. Currently, it supports the
        "uniform", "group_proportional" or "proportional" methods.

    Returns:
        numpy.ndarray: Matrix of initial probabilities
    """
    M_size, I_size = X.shape
    G_size = b.shape[1]
    if p_method == "uniform":
        p_est = np.full((G_size, I_size), 1 / I_size)
    if p_method == "proportional":
        p_est = np.array([np.sum(X, axis=0) / np.sum(X) for g in range(G_size)])
    if p_method == "group_proportional":
        p_mesas = X / np.sum(X, axis=1)[..., None]
        p_est = np.zeros((G_size, I_size))
        agg_b = np.sum(b, axis=0)
        for g in range(G_size):
            for i in range(I_size):
                p_est[g, i] = np.sum(p_mesas[:, i] * b[:, g]) / np.sum(b[:, g])
        p_est[np.isnan(p_est)] = 0
    return p_est


def compute_Tm_list(n_m, p, b_m, I_size, G_size):
    K = [[] for f in range(G_size + 1)]
    K_dict = {}
    H = [[] for f in range(G_size + 1)]
    T = [[] for f in range(G_size + 1)]
    # U = [[] for f in range(G_size+1)]
    K[-1] = np.array(combinations(I_size, 0))
    for k_ind, k in enumerate(K[-1]):
        K_dict[-1, tuple(k)] = k_ind
    T[-1] = np.full(len(K[-1]), 1)
    # U[-1] = np.full((len(K[-1]),G_size,I_size), 1)
    lgac_n = [
        sum([np.log(max(k, 1)) for k in range(n + 1)]) for n in range(np.max(b_m) + 1)
    ]
    log_p = np.log(p)
    hlg_p_gin = [
        [[n * log_p[g, i] for n in range(np.max(b_m) + 1)] for i in range(I_size)]
        for g in range(G_size)
    ]
    hb_gn = [
        [n / b_m[g] if b_m[g] > 0 else 0.0 for n in range(max(b_m) + 1)]
        for g in range(G_size)
    ]
    b_m_cum = [np.sum(b_m[: f + 1]) for f in range(G_size)]
    for f in range(G_size):
        K[f] = np.array(combinations_filtered(I_size, b_m_cum[f], n_m))
        for k_ind, k in enumerate(K[f]):
            K_dict[f, tuple(k)] = k_ind
        H[f] = np.array(combinations_filtered(I_size, b_m[f], n_m))
        T[f] = np.zeros(len(K[f]))
        # U[f] = np.zeros((len(K[f]),f+1,I_size))

        for k_ind in range(len(K[f])):
            k = K[f][k_ind]
            T_k = 0.0
            U_k = np.zeros((f + 1, I_size))
            for h_ind in range(len(H[f])):
                h = H[f][h_ind]
                if all(h <= k):
                    # k_ind_before = find_tuple(K[f-1],k-h)
                    if f == 0:
                        k_ind_before = 0
                    else:
                        k_ind_before = K_dict[f - 1, tuple(k - h)]
                        # k_ind_before = find_tuple(K[f-1],k-h)
                        # k_ind_before = get_index(k-h,n_m,b_m_cum[f-1],I_size)
                    a = np.exp(
                        lgac_n[b_m[f]]
                        + np.sum(
                            [
                                hlg_p_gin[f][i][h[i]] - lgac_n[h[i]]
                                for i in range(I_size)
                            ]
                        )
                    )
                    T_k += T[f - 1][k_ind_before] * a
                    # for i in range(I_size):
                    #     a_h_b = a * hb_gn[f][h[i]]
                    #     for g in range(f):
                    #         U_k[g,i] += U[f-1][k_ind_before,g,i]*a
                    #     if h[i] > 0:
                    #         U_k[f,i] += T[f-1][k_ind_before]*a_h_b

                T[f][k_ind] = T_k
                # U[f][k_ind] = U_k

    return T[G_size - 1][0]


def compute_ll(n, p, b):
    # sum over ln(T) for each ballot_box
    ll = 0
    for m in range(n.shape[0]):
        ll += np.log(compute_Tm_list(n[m], p, b[m], n.shape[1], b.shape[1]))

    return ll


# (g,i) compute estimate of p using EM algorithm with parameters X and b
def EM_algorithm_old(
    X,
    b,
    p_est,
    q_method,
    convergence_value=0.0001,
    max_iterations=100,
    load_bar=True,
    verbose=True,
    dict_results={},
    save_dict=False,
    dict_file=None,
):
    M_size, I_size = X.shape
    G_size = b.shape[1]
    J_mean = np.round(np.mean(np.sum(X, axis=1)), 1)
    if verbose:
        print(
            "M =",
            M_size,
            " G =",
            G_size,
            " I =",
            I_size,
            " J =",
            J_mean,
            " delta =",
            convergence_value,
        )

    if verbose:
        print("-" * 100)
        print("EM-algorithm")
        print("-" * 100)

    run_time = 0
    ## initial dict ##
    dict_results["p_est"] = p_est
    dict_results["end"] = -1
    dict_results["time"] = run_time
    dict_results["iterations"] = 0
    # save dict as pickle
    if save_dict:
        with open(dict_file, "wb") as handle:
            pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

    previous_Q = -np.inf
    for i in tqdm(range(1, max_iterations + 1), disable=not load_bar):
        start_iteration = time.time()
        q = q_method(X, p_est, b)
        p_new = compute_p(q, b)
        p_new[np.isnan(p_new)] = 0
        end_iterartion = time.time()
        # update time
        run_time += end_iterartion - start_iteration
        # update Q
        log_p_new = np.where(p_new > 0, np.log(p_new), 0)
        Q = np.sum(b * np.sum(q * log_p_new, axis=2))
        dif_Q = Q - previous_Q
        dict_results["Q"] = Q
        dict_results["dif_Q"] = dif_Q

        if verbose:
            print("-" * 50)
            print(np.round(p_new, 4))
            print("Δ: ", np.max(np.abs(p_new - p_est)))
            print("Q: ", Q)
            # print('a: ',np.sum(q * log_p_est, axis = 2))
            print("-" * 50)
        # print(np.round(p_est[1,6],5))

        # check convergence of p
        if (np.abs(p_new - p_est) < convergence_value).all():
            log_q = np.where(q > 0, np.log(q), 0)
            # changed compute conditional
            E_log_q = np.sum(b * np.sum(log_q * q, axis=2))
            dict_results["q"] = q
            dict_results["E_log_q"] = E_log_q

            dict_results["end"] = 1
            if verbose:
                print(f"Convergence took {i} iterations and {run_time} seconds.")

            # save results for convergence
            dict_results["p_est"] = p_new
            dict_results["time"] = run_time
            dict_results["iterations"] = i
            if save_dict:
                with open(dict_file, "wb") as handle:
                    pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
            return p_new, i, run_time

        # check if the expected Likelihood is not increasing
        if dif_Q < 0:
            dict_results["end"] = 2
            print("early exit! logLL decreased at iteration {}".format(i))
            if verbose:
                print(f"Q diminished; took {i} iterations and {run_time} seconds.")

            # save results for convergence
            dict_results["p_est"] = p_est  # previous one was better
            dict_results["time"] = run_time
            dict_results["iterations"] = i
            if save_dict:
                with open(dict_file, "wb") as handle:
                    pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
            return p_est, i, run_time
        previous_Q = Q

        # save results for iteration
        dict_results["p_est"] = p_new
        dict_results["end"] = 0
        dict_results["time"] = run_time
        dict_results["iterations"] = i
        if save_dict:
            with open(dict_file, "wb") as handle:
                pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

        p_est = p_new.copy()
    if verbose:
        print(f"Did not reach convergence after {i} iterations and {run_time}. \n")
    return p_new, i, run_time


## VERSION NUEVA
def EM_algorithm(
    X,
    b,
    p_est,
    q_method,
    convergence_value=0.001,
    max_iterations=100,
    load_bar=True,
    verbose=True,
    dict_results={},
    save_dict=False,
    dict_file=None,
):
    """
    Implements the whole EM algorithm.

    Parameters:
        X (numpy.ndarray): Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
        b (numpy.ndarray): Matrix of dimension (bxg) that stores the amount of votes from demographic group "g".
        p_est (numpy.ndarray): Matrix of initial probabilities.
        q_method (string): Method for estimating the probability that a voter of group "g" in ballot box "b" voted for candidate "c" conditional on the observed result. Currently, it supports
        the "", "", "" and "" methods.
        convergence_value (float, optional): The epsilon value of convergence.
            default value: 0.001
        max_iterations (int, optional): The maximum amount of iterations.
            default value: 100
        load_bar (bool, optional): Print a progress bar of the process.
            default value: True
        verbose (bool, optional): Print indicating messages.
            default value: True
        dict_results (dict, optional): Dictionary that stores the progress of the algorithm, including the initial parameters, ending criteria, run time and amount of iterations.
            default value: {}
        save_dict (bool, optional): Save the dictionary that stores the progress of the algorithm.
            default value: False
        dict_file (str, optional): The file extension of the resulting file.
            default value: None.




    """
    logLikelys = np.empty(max_iterations)
    M_size, I_size = X.shape
    G_size = b.shape[1]
    J_mean = np.round(np.mean(np.sum(X, axis=1)), 1)
    if verbose:
        print(
            "M =",
            M_size,
            " G =",
            G_size,
            " I =",
            I_size,
            " J =",
            J_mean,
            " delta =",
            convergence_value,
        )

    if verbose:
        print("-" * 100)
        print("EM-algorithm")
        print("-" * 100)

    run_time = 0
    ## initial dict ##
    dict_results["p_est"] = p_est
    dict_results["end"] = -1
    dict_results["time"] = run_time
    dict_results["iterations"] = 0
    # save dict as pickle
    if save_dict:
        with open(dict_file, "wb") as handle:
            pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

    previous_Q = -np.inf
    previous_ll = -np.inf
    dict_results["end"] = 0  # not converged yet until the time/iteration limit
    start_iteration = (
        perf_counter()
    )  # This iteration starts before the loop NOW, to compare it against C.
    for i in tqdm(range(max_iterations), disable=not load_bar):
        # start_iteration = time.time()
        q = q_method(X, p_est, b)
        p_new = compute_p(q, b)
        p_new[np.isnan(p_new)] = 0

        # check convergence of p

        """
        if (np.abs(p_new - p_est) < convergence_value).all():
            end_time = perf_counter()
            run_time = end_time - start_iteration
            # run_time += perf_counter() - start_iteration
            dict_results["end"] = 1
            if verbose:
                print(f"Convergence took {i} iterations and {run_time} seconds.")
"""
        # update time
        # run_time += perf_counter() - start_iteration

        # update Q
        log_p_new = np.where(p_new > 0, np.log(p_new + epsilon), 0)
        Q = np.sum(b * np.sum(q * log_p_new, axis=2))
        dif_Q = Q - previous_Q
        dict_results["dif_Q"] = dif_Q
        dict_results["Q"] = Q

        # compute the other term of the expected log likelihood
        log_q = np.where(q > 0, np.log(q + epsilon), 0)
        E_log_q = np.sum(b * np.sum(log_q * q, axis=2))
        dict_results["q"] = q
        dict_results["E_log_q"] = E_log_q

        # save expected log likelihood
        # dict_results["ll"] = Q - dict_results["E_log_q"]
        dict_results["ll"] = compute_ll(X, p_new, b)  # new version
        logLikelys[i] = dict_results["ll"]

        if verbose:
            print("-" * 50)
            print("iteration: ", i)
            print(np.round(p_new, 4))
            print("Δ: ", np.max(np.abs(p_new - p_est)))
            print("Q: ", Q)
            print("ll: ", dict_results["ll"])
            # print('a: ',np.sum(q * log_p_est, axis = 2))
            print("-" * 50)
        # print(np.round(p_est[1,6],5))

        # check if the expected Likelihood is not increasing
        """
        if previous_ll - dict_results["ll"] > 0:
            p_new = p_est.copy()
            dict_results["end"] = 2
            if verbose:
                print(f"ll decreased; took {i} iterations and {run_time} seconds.")
        """
        previous_ll = dict_results["ll"].copy()

        # save results for iteration
        dict_results["p_est"] = p_new

        dict_results["time"] = run_time
        dict_results["iterations"] = i
        if save_dict:
            with open(dict_file, "wb") as handle:
                pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # update p
        p_est = p_new.copy()

        if dict_results["end"] > 0:
            break

    if save_dict:
        with open(dict_file, "wb") as handle:
            pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if dict_results["end"] == 0:
        if verbose:
            print(f"Did not reach convergence after {i} iterations and {run_time}. \n")

    # fix if one group doesnt have voters
    agg_b = np.sum(b, axis=0)
    p_new[agg_b == 0, :] = np.sum(X, axis=0) / np.sum(X)

    return p_new, i, run_time, logLikelys


# (g,i) compute estimate of p using EM algorithm with parameters X and b
