import numpy as np
from time import perf_counter
from tqdm import tqdm
import pickle

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
    I_size = X.shape
    G_size = b.shape[1]

    if p_method == "uniform":
        p_est = np.full((G_size, I_size), 1 / I_size)

    if p_method == "proportional":
        p_est = np.array([np.sum(X, axis=0) / np.sum(X) for g in range(G_size)])

    if p_method == "group_proportional":
        p_mesas = X / np.sum(X, axis=1)[..., None]  # Sum over the ballot axis
        p_est = np.zeros((G_size, I_size))  # GxI

        for g in range(G_size):
            for i in range(I_size):
                p_est[g, i] = np.sum(p_mesas[:, i] * b[:, g]) / np.sum(
                    b[:, g]
                )  # Maybe an improvement would be to have the sum over the ballot instead of universe

        p_est[np.isnan(p_est)] = (
            0  # This is to avoid border cases where there's division by cero.
        )

    return p_est


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
    dict_results["end"] = -1  # -1 is used when the stop is not made by convergence
    dict_results["time"] = run_time
    dict_results["iterations"] = 0

    # save dict as pickle
    if save_dict:
        with open(dict_file, "wb") as handle:
            pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

    previous_Q = -np.inf
    previous_ll = -np.inf
    dict_results["end"] = 0  # 0 if it hasn't converged until the time/iteration limit.

    for i in tqdm(range(1, max_iterations + 1), disable=not load_bar):
        start_iteration = perf_counter()
        q = q_method(X, p_est, b)  # Calculates q with method of choice
        p_new = compute_p(q, b)
        p_new[np.isnan(p_new)] = 0

        # Check convergence of p
        if (np.abs(p_new - p_est) < convergence_value).all():
            dict_results["end"] = 1
            if verbose:
                print(f"Convergence took {i} iterations and {run_time} seconds.")

        # Update time
        run_time += perf_counter() - start_iteration

        # Update Q; this is for evaluation indicators
        log_p_new = np.where(p_new > 0, np.log(p_new), 0)
        Q = np.sum(b * np.sum(q * log_p_new, axis=2))
        dif_Q = Q - previous_Q
        dict_results["dif_Q"] = dif_Q
        dict_results["Q"] = Q

        # Compute the other term of the expected log likelihood; for indicators too
        log_q = np.where(q > 0, np.log(q), 0)
        E_log_q = np.sum(b * np.sum(log_q * q, axis=2))
        dict_results["q"] = q
        dict_results["E_log_q"] = E_log_q

        # Save expected log likelihood
        dict_results["ll"] = Q - dict_results["E_log_q"]

        if verbose:
            print("-" * 50)
            print("iteration: ", i)
            print(np.round(p_new, 4))
            print("Δ: ", np.max(np.abs(p_new - p_est)))
            print("Q: ", Q)
            print("ll: ", dict_results["ll"])
            print("-" * 50)

        # Check if the expected likelihood is not increasing:
        # Border case, seemingly an decimal approximation error, where the LL starts decreasing and
        # iterating through two values.

        if previous_ll - dict_results["ll"] > 0:
            p_new = p_est.copy()
            dict_results["end"] = 2
            if verbose:
                print(f"ll decreased; took {i} iterations and {run_time} seconds.")
        previous_ll = dict_results["ll"].copy()

        # Save results for iteration
        dict_results["p_est"] = p_new
        dict_results["time"] = run_time
        dict_results["iterations"] = i
        if save_dict:
            with open(dict_file, "wb") as handle:
                pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # Update p for the next iteration
        p_est = p_new.copy()  # TODO: Make it a pointer in C

        if dict_results["end"] > 0:
            break

    if save_dict:
        with open(dict_file, "wb") as handle:
            pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if dict_results["end"] == 0:
        if verbose:
            print(f"Did not reach convergence after {i} iterations and {run_time}. \n")

    # Fix if one group doesnt have voters; border case.
    agg_b = np.sum(b, axis=0)
    p_new[agg_b == 0, :] = np.sum(X, axis=0) / np.sum(X)

    return p_new, i, run_time
