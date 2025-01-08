import numpy as np
import time
from tqdm import tqdm
from EM_algorithm import EM_algorithm, get_p_est
from helper_functions import readMatrixFromC
from time import perf_counter

epsilon = 1e-10  # Small stabilizing value


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
    print("Running the EM algorithm with the multinomial method")
    matrixList = readMatrixFromC(
        "project/src/matricesTest3.bin"
    )  # First one is X, second one is W/b

    start_iteration2 = (
        perf_counter()
    )  # This iteration starts before the loop NOW, to compare it against C.
    answer = EM_mult(matrixList[0].T, matrixList[1], max_iterations=10000, verbose=True)
    end_time = perf_counter()
    run_time = end_time - start_iteration2
    print("El tiempo de ejecución es {}".format(run_time))
    print(answer[0], answer[1], answer[2])
