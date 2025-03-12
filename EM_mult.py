import numpy as np
import time
from tqdm import tqdm
from EM_algorithm import EM_algorithm, get_p_est
from helper_functions import readMatrixFromC
from time import perf_counter
import threading

epsilon = 1e-10  # Small stabilizing value


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


def compute_q_multinomial(n, p, b):
    p_cond = (b @ p)[:, None, :] - p[None, ...]
    p_cond = np.where(p_cond == 0, epsilon, p_cond)  # Replace zeros with epsilon
    q = p[None, ...] * np.expand_dims(n, axis=1) / p_cond
    q[np.isnan(q)] = 0
    sum_q = np.sum(q, axis=2)[:, :, None]
    sum_q = np.where(sum_q == 0, epsilon, sum_q)  # Replace zeros with epsilon
    q = q / sum_q
    # q = q / np.sum(q, axis=2)[:, :, None]
    q[np.isnan(q)] = 0
    return q


# def compute_q_multinomial(n,p,b):
#     p_cond = (b @ p)[:,None,:] - p[None,...]
#     q = np.expand_dims(n, axis=1) / p_cond
#     # print('r')
#     # print(p_cond)
#     # print('q')
#     # print(q/np.sum(q, axis=2)[:,:,None])
#     # print('p')
#     # print(p[None,...])
#     # exit()
#     q[np.isnan(q)] = 0
#     q = p[None,...]*q/np.sum(q, axis=2)[:,:,None]
#     #q[np.isnan(q)] = 0
#     return q


def compute_p(q, b):
    num = np.sum(np.multiply(q, b[..., None]), axis=0)
    dem = np.sum(b, axis=0)[..., None]
    return num / dem


# (g,i) compute estimate of p using EM algorithm with parameters X and b
def EM_mult(
    X,
    b,
    p_est=None,
    convergence_value=0.0001,
    max_iterations=100,
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
        compute_q_multinomial,
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
plt.title("Multinomial: Average Log-Likelihood per Iteration (20 seeds) (C2G2)")
# plt.ylim(-45, 150)
plt.legend()
plt.grid(True, which="both")
plt.show()
