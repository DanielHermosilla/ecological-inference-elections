from time import perf_counter
import numpy as np
from scipy.stats import multivariate_normal, multinomial
from helper_functions import *
import threading
import matplotlib.pyplot as plt
import json

# from model_elections import compute_qm_list
from multiprocessing import Pool

# from helper_functions import combinations
import time
from tqdm import tqdm
from EM_algorithm import EM_algorithm, get_p_est


def process_instance(s, loglikelihood_list):
    """Function to process each instance in a separate thread."""
    with open(f"instances/J100_M50_G2_I2_L50_seed{s}.json") as f:
        instance = json.load(f)

    X = np.array(instance["X"])
    b = np.array(instance["W"])
    p = np.array(instance["p"])

    start_time = perf_counter()

    answer = EM_mvn_pdf(
        X,
        b,
        convergence_value=0.001,
        p_method="group_proportional",
        max_iterations=1000,
        verbose=False,
        load_bar=True,
    )

    end_time = perf_counter()
    run_time = end_time - start_time

    loglikelihood_list[s - 1] = answer[3]  # Store result in correct position


def compute_q_mvn_pdf(n, p, b, parallel=False):
    M_size, G_size, I_size = b.shape[0], b.shape[1], n.shape[1]
    q = np.zeros(shape=(M_size, G_size, I_size))
    p_trunc = p[:, :-1]
    n_trunc = n[:, :-1]
    diag_p = [np.diag(p_g) for p_g in p_trunc]
    p_g_squared = np.einsum("ij,ik->ijk", p_trunc, p_trunc)
    # if M_size >= 100:
    #     parallel = True
    # if parallel:
    #     q_args = [(n[m], p, b[m], I_size, G_size, use_pdf) for m in range(M_size)]
    #     with Pool() as pool:
    #         q = pool.starmap(compute_qm_mvn, q_args)
    #     return np.array(q)
    for m in range(M_size):
        q_m = compute_qm_mvn_pdf(
            n_trunc[m], p, p_trunc, b[m], I_size, G_size, diag_p, p_g_squared
        )
        q[m] = q_m
    return q


def compute_qm_mvn_pdf(n_m_trunc, p, p_trunc, b_m, I_size, G_size, diag_p, p_g_squared):
    """
    Notation:
    - b_m = The group votes of ballot "m". (1xGROUP_SIZE)
    - p_trunc = The probabilities for the group "g" and candidate "c", without considering the last candidate.

    """
    # It will omit the last candidate since it's redundant.
    mu = (
        b_m @ p_trunc
    )  # (1,G_size) @ (G_size, I_size-1) = I_size - 1 # Calculates the mean, not considering the last candidate.
    # Covariate adjusted matrix, which corresponds as;
    # $$\Sigma^{g}_{b}=diag(\mu_{b}^{g})-p^{T}_{\text{trunc}}\cdot diag(b_m)\cdot p_{\text{trunc}}$$
    cov = np.diag(mu) - p_trunc.T @ np.diag(b_m) @ p_trunc  # (I_size-1, I_size-1)

    # Corrects the covariance, centers it. The same for the average.
    covs_U = cov - diag_p + p_g_squared  # (G_size, I_size-1, I_size-1)
    # print("The first covariance matrix is:")
    # print(covs_U)
    mus_U = mu - p_trunc  # (G_size, I_size-1)
    # print("The first mu matrix is:")
    # print(mus_U)

    vals_U, vecs_U = np.linalg.eigh(covs_U)

    mahas = np.zeros(shape=(G_size, I_size))

    inverses, invs_devs, mahas[:, I_size - 1] = get_maha_mult_v2(
        n_m_trunc, mus_U, vals_U, vecs_U
    )

    diag_inverses = np.diagonal(inverses, axis1=1, axis2=2)
    # print("The multiplication from mahanobis is")
    # print(invs_devs)
    mahas[:, :-1] = mahas[:, I_size - 1][..., None] - 2 * invs_devs + diag_inverses
    # print("The mahanalobis values are:")
    # print(mahas)

    q_m = np.exp(-0.5 * mahas) * p  # agregar mahalanobis de T
    q_m = q_m / np.sum(q_m, axis=1)[:, None]

    return q_m


# https://gregorygundersen.com/blog/2019/10/30/scipy-multivariate/
# https://gregorygundersen.com/blog/2020/12/12/group-multivariate-normal-pdf/
def get_maha_mult_v2(x, means, vals, vecs):
    # Invert the eigenvalues.
    valsinvs = 1.0 / vals

    devs = x - means  # G x I - 1
    # print("The diferences are")
    # print(devs)

    valsinv_diag = [np.diag(v) for v in valsinvs]

    inverses = vecs @ valsinv_diag @ vecs.swapaxes(1, 2)  # G x I-1 x I-1
    # print("The inverse matrix is")
    # print(inverses)
    invs_devs = np.einsum(
        "gi,gij->gj", devs, inverses
    )  # G x I-1 (for each g you get the left hand side of mahalanobis)

    mahas = np.einsum(
        "gi,gi->g", devs, invs_devs
    )  # G (for each g you get the mahalanobis distance)

    return inverses, invs_devs, mahas


def compute_p(q, b):
    num = np.sum(np.multiply(q, b[..., None]), axis=0)
    dem = np.sum(b, axis=0)[..., None]
    return num / dem


def EM_mvn_pdf(
    X,
    b,
    p_est=None,
    convergence_value=0.0001,
    max_iterations=10000,
    p_method="group_proportional",
    load_bar=True,
    verbose=True,
    dict_results={},
    save_dict=False,
    dict_file=None,
):
    if p_est is None:
        p_est = get_p_est(X, b, p_method)
    return EM_algorithm(
        X,
        b,
        p_est,
        compute_q_mvn_pdf,
        convergence_value,
        max_iterations,
        load_bar,
        verbose,
        dict_results,
        save_dict,
        dict_file,
    )


if __name__ == "__main__":
    num_threads = 20  # Number of JSON files to process
    loglikelihood_list = [None] * num_threads  # Pre-allocate storage for results
    threads = []

    start_global = perf_counter()  # Global timer

    # Create and start threads
    for s in range(1, num_threads + 1):
        t = threading.Thread(target=process_instance, args=(s, loglikelihood_list))
        threads.append(t)
        t.start()

    # Wait for all threads to finish
    for t in threads:
        t.join()

    end_global = perf_counter()
    print(f"Total execution time: {end_global - start_global:.2f} seconds")
    decreasing_info = {}

    for i, loglikelihood in enumerate(loglikelihood_list):
        decreasing_steps = np.where(np.diff(loglikelihood) < 0)[
            0
        ]  # Encuentra índices donde decrece
        if len(decreasing_steps) > 0:
            decreasing_info[i] = decreasing_steps  # Guardamos en un diccionario

    # Mostrar resultados
    for thread_idx, indices in decreasing_info.items():
        print(f"Thread {thread_idx + 1} decrece en las posiciones: {indices}")

    # Convert loglikelihood results to NumPy array
    loglikelihood_matrix = np.array(loglikelihood_list)

    # Compute mean log-likelihood per iteration
    loglikelihood_mean = np.mean(loglikelihood_matrix, axis=0)
    iterations = np.arange(1, 1001)  # Avoid zero, since log(0) is undefined
iterations2 = np.arange(1, 1000)  # Avoid zero, since log(0) is undefined

# Plot
plt.figure(figsize=(8, 5))
plt.plot(
    iterations,
    loglikelihood_mean,
    color="blue",
    alpha=0.5,
)
plt.plot(iterations, loglikelihood_mean, "ro", color="red", alpha=0.7, markersize=1)
for i in range(num_threads):
    plt.plot(iterations, loglikelihood_list[i], label=f"(S{i+1})", alpha=0.2)
    plt.plot(iterations2, np.diff(loglikelihood_list[i]))
plt.xlabel("Iteration")
plt.ylabel("Log-Likelihood")
plt.xscale("log")
plt.title("PDF: Average Log-Likelihood per Iteration (20 seeds) (C2G2)")
# plt.ylim(-45, 150)
plt.legend()
plt.grid(True, which="both")
plt.show()
