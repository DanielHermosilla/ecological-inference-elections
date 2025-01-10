from time import perf_counter
import numpy as np
from scipy.stats import multivariate_normal, multinomial
from helper_functions import *

# from model_elections import compute_qm_list
from multiprocessing import Pool

# from helper_functions import combinations
import time
from tqdm import tqdm
from EM_algorithm import EM_algorithm, get_p_est


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
    print("Running the multivariate PDF algorithm")
    matrixList = readMatrixFromC(
        "project/src/matricesTest3.bin"
    )  # First one is X, second one is W/b

    start_iteration2 = (
        perf_counter()
    )  # This iteration starts before the loop NOW, to compare it against C.
    answer = EM_mvn_pdf(
        matrixList[0].T, matrixList[1], max_iterations=200, verbose=False
    )
    end_time = perf_counter()
    run_time = end_time - start_iteration2
    print("El tiempo de ejecución es {}".format(run_time))
    print(answer[0], answer[1], answer[2])
