import numpy as np
from scipy.stats import multivariate_normal, multinomial

# from model_elections import compute_qm_list
from multiprocessing import Pool

# from helper_functions import combinations
import time
from tqdm import tqdm

# import rpy2.robjects as robjects
from EM_algorithm import EM_algorithm, get_p_est

# from EM_full import EM_full
from scipy.stats import norm

# r = robjects.r
# r_code = """
# library(mvtnorm)
# """
# r(r_code)

# from EM_algorithm import EM_algorithm, get_p_est
# from EM_full import EM_full, compute_q_list
# from EM_mvn_pdf import EM_mvn_pdf, compute_q_mvn_pdf
# from EM_mult import EM_mult, compute_q_multinomial


def compute_q_mvn_cdf(n, p, b, parallel=False, R=False):
    """
    Computes the matrix with the probabilities that a voter of group "g" in ballot box "b" voted for candidate "c" (referred as "q" on the paper).

    Parameters:
        n (numpy.ndarray): Matrix with the amount of votes for each candidate (referred as "X" on the paper)
        p (numpy.ndarray): Matrix with the prior M-step probabilities.
        b (numpy.ndarray): Matrix with the amount of votes per demographic group.

    Returns:
        numpy.ndarray: Matrix with all the values of "q".
    """
    M_size, G_size, I_size = b.shape[0], b.shape[1], n.shape[1]
    q = np.zeros(shape=(M_size, G_size, I_size))
    p_trunc = p[
        :, :-1
    ]  # For simplicity, removes the last multivariate normal (maybe include another method for dropping redundant MVN)
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
        q_m = compute_qm_mvn_cdf(
            n[m], n_trunc[m], p, p_trunc, b[m], I_size, G_size, diag_p, p_g_squared
        )
        q[m] = q_m
    return q


def compute_qm_mvn_cdf(
    n_m, n_m_trunc, p, p_trunc, b_m, I_size, G_size, diag_p, p_g_squared
):

    mu = b_m @ p_trunc  # (1,G_size) @ (G_size, I_size-1) = I_size - 1
    cov = np.diag(mu) - p_trunc.T @ np.diag(b_m) @ p_trunc  # (I_size-1, I_size-1)

    covs_U = cov - diag_p + p_g_squared  # (G_size, I_size-1, I_size-1)
    mus_U = mu - p_trunc  # (G_size, I_size-1)

    if I_size > 2:
        Chols_U = np.linalg.cholesky(covs_U)

    qm = np.zeros(shape=(G_size, I_size))
    for i in range(I_size):
        n_i = n_m_trunc.copy()
        if i < I_size - 1:
            n_i[i] -= 1
        for g in range(G_size):
            if I_size == 2:
                qm[g, i] = norm.cdf(n_i + 0.5, mus_U[g], np.sqrt(covs_U[g])) - norm.cdf(
                    n_i - 0.5, mus_U[g], np.sqrt(covs_U[g])
                )
            else:
                n_i_center = n_i - mus_U[g]
                qm[g, i] = MonteCarlo_cdf_matrix(
                    Chols_U[g], n_i_center - 0.5, n_i_center + 0.5, I_size - 1, 0.0001
                )
    if np.all(qm == 0):
        qm = np.ones(shape=(G_size, I_size))
    # print(qm)
    # print('----')
    qm = qm * p  # agregar término p
    qm = qm / np.sum(qm, axis=1)[:, None]  # normalize
    # check if nan
    # if np.isnan(qm).any():
    #     exit()
    return qm


# def compute_qm_mvn_cdf_R(n_m, n_m_trunc, p, p_trunc, b_m, I_size, G_size, diag_p, p_g_squared, s = 100000):
#     mu = b_m @ p_trunc   # (1,G_size) @ (G_size, I_size-1) = I_size - 1
#     cov = np.diag(mu) - p_trunc.T @ np.diag(b_m) @ p_trunc # (I_size-1, I_size-1)

#     covs_U = cov - diag_p + p_g_squared # (G_size, I_size-1, I_size-1)
#     mus_U = mu - p_trunc # (G_size, I_size-1)

#     qm = np.zeros(shape=(G_size, I_size))
#     for i in range(I_size):
#         n_i = n_m_trunc.copy()
#         if i < I_size - 1:
#             n_i[i] -= 1
#         for g in range(G_size):
#             mean_vector = robjects.FloatVector(mus_U[g])
#             cov_matrix = robjects.r.matrix(covs_U[g], nrow=I_size-1, ncol=I_size-1)
#             lower_bound = robjects.FloatVector(n_i - 0.5)
#             upper_bound = robjects.FloatVector(n_i + 0.5)
#             qm_gi = r.mvtnorm.pmvnorm(lower_bound, upper_bound, mean_vector, cov_matrix)
#             print(qm_gi)
#             qm[g,i] = qm_gi
#             # X = np.random.multivariate_normal(mus_U[g], covs_U[g], s)
#             # inside_cube = (X >= n_i - 0.5) * (X <= n_i + 0.5)
#             # fav_cases = np.all(inside_cube, axis=1)
#             # qm[g,i] = np.sum(fav_cases)/s

#     qm = qm * p # agregar término p
#     qm = qm/np.sum(qm, axis=1)[:,None] # normalize

#     return qm


def MonteCarlo_cdf(Chol, a, b, mvn_dim, epsilon, alpha, Nmax):
    intsum = 0
    N = 1
    varsum = 0
    d = np.zeros(mvn_dim)
    e = np.zeros(mvn_dim)
    f = np.zeros(mvn_dim)
    y = np.zeros(mvn_dim)
    d[0] = norm.cdf(a[0] / Chol[0, 0])
    e[0] = norm.cdf(b[0] / Chol[0, 0])
    f[0] = e[0] - d[0]
    error = 1
    while error > epsilon and N < Nmax:
        w = np.random.uniform(0, 1, mvn_dim - 1)
        for i in range(1, mvn_dim):
            y[i - 1] = norm.ppf(d[i - 1] + w[i - 1] * (e[i - 1] - d[i - 1]))
            Chol_cum = np.dot(Chol[i, :i], y[:i])
            d[i] = norm.cdf((a[i] - Chol_cum) / Chol[i, i])
            e[i] = norm.cdf((b[i] - Chol_cum) / Chol[i, i])
            f[i] = (e[i] - d[i]) * f[i - 1]
        N = N + 1
        # print(f[mvn_dim-1])
        varsum = varsum + f[mvn_dim - 1] ** 2
        intsum = intsum + f[mvn_dim - 1]
        error = alpha * np.sqrt((varsum / N - (intsum / N) ** 2) / N)
    return intsum / N


def MonteCarlo_cdf(Chol, a, b, mvn_dim, epsilon, Nmax):
    intsum = 0
    N = 1
    varsum = 0
    d = np.zeros(mvn_dim)
    e = np.zeros(mvn_dim)
    f = np.zeros(mvn_dim)
    y = np.zeros(mvn_dim)
    d[0] = norm.cdf(a[0] / Chol[0, 0])
    e[0] = norm.cdf(b[0] / Chol[0, 0])
    f[0] = e[0] - d[0]
    error = 1
    while error > epsilon and N < Nmax:
        w = np.random.uniform(0, 1, mvn_dim - 1)
        for i in range(1, mvn_dim):
            y[i - 1] = norm.ppf(d[i - 1] + w[i - 1] * (e[i - 1] - d[i - 1]))
            Chol_cum = np.dot(Chol[i, :i], y[:i])
            d[i] = norm.cdf((a[i] - Chol_cum) / Chol[i, i])
            e[i] = norm.cdf((b[i] - Chol_cum) / Chol[i, i])
            f[i] = (e[i] - d[i]) * f[i - 1]
        N = N + 1
        # print(f[mvn_dim-1])
        varsum = varsum + f[mvn_dim - 1] ** 2
        intsum = intsum + f[mvn_dim - 1]
        error = np.sqrt((varsum / N - (intsum / N) ** 2) / N)
    return intsum / N


def MonteCarlo_cdf_matrix(Chol, a, b, mvn_dim, epsilon, min_order=2, max_order=5):
    intsum = 0
    varsum = 0
    total_sim = 0
    N = [10**i for i in range(min_order, max_order + 1)]
    d_0 = norm.cdf(a[0] / Chol[0, 0])
    e_0 = norm.cdf(b[0] / Chol[0, 0])
    f_0 = e_0 - d_0
    for n in N:
        d = np.zeros((mvn_dim, n))
        e = np.zeros((mvn_dim, n))
        f = np.zeros((mvn_dim, n))
        y = np.zeros((mvn_dim, n))
        d[0, :] = d_0
        e[0, :] = e_0
        f[0, :] = f_0
        w = np.random.uniform(0, 1, (mvn_dim - 1, n))
        for i in range(1, mvn_dim):
            y[i - 1, :] = norm.ppf(
                d[i - 1, :] + w[i - 1, :] * (e[i - 1, :] - d[i - 1, :])
            )
            Chol_cum = np.einsum("j,js->s", Chol[i, :i], y[:i, :])
            d[i, :] = norm.cdf((a[i] - Chol_cum) / Chol[i, i])
            e[i, :] = norm.cdf((b[i] - Chol_cum) / Chol[i, i])
            f[i, :] = (e[i, :] - d[i, :]) * f[i - 1, :]
        total_sim += n
        varsum = varsum + np.sum(f[mvn_dim - 1, :] ** 2)
        intsum = intsum + np.sum(f[mvn_dim - 1, :])
        error = np.sqrt((varsum / total_sim - (intsum / total_sim) ** 2) / total_sim)
        if error < epsilon:
            break

    return intsum / total_sim

    # d[0] = norm.cdf(a[0]/Chol[0,0])
    # e[0] = norm.cdf(b[0]/Chol[0,0])
    # f[0] = e[0] - d[0]
    # error = 1
    # while error > epsilon and N < Nmax:
    #     w = np.random.uniform(0,1,mvn_dim-1)
    #     for i in range(1,mvn_dim):
    #         y[i-1] = norm.ppf(d[i-1] + w[i-1]*(e[i-1]-d[i-1]))
    #         Chol_cum = np.dot(Chol[i,:i],y[:i])
    #         d[i] = norm.cdf((a[i]-Chol_cum)/Chol[i,i])
    #         e[i] = norm.cdf((b[i]-Chol_cum)/Chol[i,i])
    #         f[i] = (e[i] - d[i])*f[i-1]
    #     N = N + 1
    #     # print(f[mvn_dim-1])
    #     varsum = varsum + f[mvn_dim-1]**2
    #     intsum = intsum + f[mvn_dim-1]
    #     error = alpha*np.sqrt((varsum/N - (intsum/N)**2)/N)
    # return intsum/N


# (g,i) compute estimate of p using EM algorithm with parameters X and b
def EM_mvn_cdf(
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
        compute_q_mvn_cdf,
        convergence_value,
        max_iterations,
        load_bar,
        verbose,
        dict_results=dict_results,
        save_dict=save_dict,
        dict_file=dict_file,
    )


if __name__ == "__main__":
    # load instance json
    import json

    with open("instances\J200_M50_G4_I2_L50_seed12.json") as f:
        instance = json.load(f)
    X = np.array(instance["n"])
    b = np.array(instance["b"])
    p = np.array(instance["p"])
    print("X ", X[-1])
    print("")
    print("b ", b[-1])
    print("")
    # compute_q_mvn_cdf(X[[0]], p, b[[0]])
    em = EM_mvn_cdf(
        X,
        b,
        convergence_value=0.001,
        p_method="uniform",
        max_iterations=100,
        verbose=True,
        load_bar=False,
    )

    # p_full = EM_full(X, b, convergence_value = 0.001, max_iterations = 100, verbose = False, load_bar = False)[0]
    # print abs error
    # print(np.sum(np.abs(p-p_est)))
    # print(np.sum(np.abs(p-p_full)))
    # print(EM_full(X,b, verbose=False, load_bar=False, p_method='group_proportional'))
    # print(EM_mvn_pdf(X,b, verbose=False, load_bar=False, p_method='group_proportional'))
    # print(EM_mult(X,b, verbose=False, load_bar=False, p_method='group_proportional'))

