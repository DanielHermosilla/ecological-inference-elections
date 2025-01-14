import numpy as np
from multiprocessing import Pool

# from poisson_multinomial import compute_q_mvn, compute_q_mvn_cdf, compute_p_value_mvn, compute_q_multinomial  # UNCOMMENT IF EXECUTING THIS FILE
import time
from time import perf_counter
from helper_functions import *
from tqdm import tqdm
from EM_algorithm import EM_algorithm, get_p_est


# @profile
def compute_qm_list(n_m, p, b_m, I_size, G_size):
    K = [[] for f in range(G_size + 1)]
    K_dict = {}
    H = [[] for f in range(G_size + 1)]
    T = [[] for f in range(G_size + 1)]
    U = [[] for f in range(G_size + 1)]
    K[-1] = np.array(combinations(I_size, 0))
    for k_ind, k in enumerate(K[-1]):
        K_dict[-1, tuple(k)] = k_ind
    T[-1] = np.full(len(K[-1]), 1)
    U[-1] = np.full((len(K[-1]), G_size, I_size), 1)
    lgac_n = [
        sum([np.log(max(k, 1)) for k in range(n + 1)])
        for n in range(int(np.max(b_m)) + 1)
    ]
    log_p = np.log(p)
    hlg_p_gin = [
        [[n * log_p[g, i] for n in range(int(np.max(b_m)) + 1)] for i in range(I_size)]
        for g in range(G_size)
    ]
    hb_gn = [
        [n / b_m[g] if b_m[g] > 0 else 0.0 for n in range(int(max(b_m)) + 1)]
        for g in range(G_size)
    ]
    b_m_cum = [np.sum(b_m[: f + 1]) for f in range(G_size)]
    for f in range(G_size):
        K[f] = np.array(combinations_filtered(I_size, int(b_m_cum[f]), n_m))
        for k_ind, k in enumerate(K[f]):
            K_dict[f, tuple(k)] = k_ind
        H[f] = np.array(combinations_filtered(I_size, int(b_m[f]), n_m))
        T[f] = np.zeros(len(K[f]))
        U[f] = np.zeros((len(K[f]), f + 1, I_size))

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
                    for i in range(I_size):
                        a_h_b = a * hb_gn[f][h[i]]
                        for g in range(f):
                            U_k[g, i] += U[f - 1][k_ind_before, g, i] * a
                        if h[i] > 0:
                            U_k[f, i] += T[f - 1][k_ind_before] * a_h_b

                T[f][k_ind] = T_k
                U[f][k_ind] = U_k

    return U[G_size - 1][0] / T[G_size - 1][0]


def compute_q_list(n, p, b, parallel=False):
    M_size, G_size, I_size = b.shape[0], b.shape[1], n.shape[1]
    q = np.zeros(shape=(M_size, G_size, I_size))
    if parallel:
        g_orders = [np.argsort(b[m]) for m in range(M_size)]
        q_args = [
            (n[m], p[g_orders[m]], b[m][g_orders[m]], I_size, G_size)
            for m in range(M_size)
        ]
        with Pool() as pool:
            q_ordered = pool.starmap(compute_qm_list, q_args)
        for m in range(M_size):
            q[m] = q_ordered[m][np.argsort(g_orders[m])]
    else:
        for m in range(M_size):
            g_order = np.argsort(b[m])  # re ordenar aplicando argsort denuevo
            back_order = np.argsort(g_order)
            q_ordered = compute_qm_list(n[m], p[g_order], b[m][g_order], I_size, G_size)
            q[m] = q_ordered[back_order]
    return q


def compute_p(q, b):
    num = np.sum(np.multiply(q, b[..., None]), axis=0)
    dem = np.sum(b, axis=0)[..., None]
    print("Current p is")
    print(num / dem)
    return num / dem


# (g,i) compute estimate of p using EM algorithm with parameters X and b
# def EM_full(X, b, convergence_value = 0.0001, max_iterations = 100,
#                  p_method = 'proportional', load_bar = True, verbose = True,
#                  dict_results = {}, save_dict = False, dict_file = None):
#     p_est = get_p_est(X, b, p_method)
#     return EM_algorithm(X, b, p_est, compute_q_list, convergence_value, max_iterations, load_bar, verbose,
#                         dict_results = dict_results, save_dict = save_dict, dict_file = dict_file)


# ADD p_est as None, now u can give it as a parameter, if None it will be computed with p_method
def EM_full(
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
        compute_q_list,
        convergence_value,
        max_iterations,
        load_bar,
        verbose,
        dict_results=dict_results,
        save_dict=save_dict,
        dict_file=dict_file,
    )


if __name__ == "__main__":
    print("Running the EM algorithm with the exact method")
    matrixList = readMatrixFromC(
        "project/src/matricesTest5.bin"
    )  # First one is X, second one is W/b

    start_iteration2 = (
        perf_counter()
    )  # This iteration starts before the loop NOW, to compare it against C.
    answer = EM_full(
        matrixList[0].T.astype(int),
        matrixList[1].astype(int),
        max_iterations=1,
        verbose=False,
    )
    end_time = perf_counter()
    run_time = end_time - start_iteration2
    print("El tiempo de ejecución es {}".format(run_time))
    print(answer[0], answer[1], answer[2])
