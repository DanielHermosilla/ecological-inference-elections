import numpy as np
import pickle
from tqdm import tqdm
import scipy.stats as stats
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor
import time


def compute_p_value_mult(n, p, b, S=100000, load_bar=False, seed=None):
    if seed is not None:
        np.random.seed(seed)
    M_size, G_size, I_size = b.shape[0], b.shape[1], n.shape[1]
    p_value_list = []
    for m in tqdm(range(M_size), disable=not load_bar):
        r_m = (p.T @ b[m]) / np.sum(b[m])
        p_value_list.append(compute_p_value_m_mult(n[m], r_m, S, seed=seed))
    return np.array(p_value_list)


def compute_p_value_m_mult(n, r, S, seed=None):
    if seed is not None:
        np.random.seed(seed)
    J = sum(n)
    lgac_n = np.array(
        [sum([np.log(max(k, 1)) for k in range(j + 1)]) for j in range(J + 1)]
    )
    log_p = np.log(r)
    x_samples = np.random.multinomial(J, r, size=S)
    beta_S = np.sum(x_samples * log_p, axis=1) - np.sum(lgac_n[x_samples], axis=1)
    beta_n = np.sum(n * log_p) - np.sum(lgac_n[n])
    less_p = np.sum(beta_S <= beta_n)
    return less_p / S


def p_val_threshold_n(n, mu, alpha):
    cum_prob = 0
    for z in range(n + 1):
        cum_prob += stats.binom.pmf(z, n, mu)
        if cum_prob >= 1 - alpha:
            return z
    return n  # In case the loop finishes without reaching threshold


def compute_thresholds(S_min, S_max, mu_power, alpha_power):
    thresholds = {}
    mu = 10 ** (-mu_power)
    alpha = 10 ** (-alpha_power)
    for s in range(S_min, S_max + 1):
        n = int(10**s)
        thresholds[s] = p_val_threshold_n(n, mu, alpha)
    return thresholds


def compute_p_value_m_mult_threshold(x, r, S_min, S_max, thresholds, seed=None):
    if seed is not None:
        np.random.seed(seed)
    J = sum(x)
    log_p = np.log(r)
    lgac_n = np.array(
        [sum([np.log(max(k, 1)) for k in range(j + 1)]) for j in range(J + 1)]
    )

    for s in range(S_min, S_max + 1):
        n = int(10**s)
        threshold = thresholds[s]
        x_samples = np.random.multinomial(J, r, size=n)
        beta_S = np.sum(x_samples * log_p, axis=1) - np.sum(lgac_n[x_samples], axis=1)
        beta_n = np.sum(x * log_p) - np.sum(lgac_n[x])
        less_p = np.sum(beta_S <= beta_n)
        if less_p >= threshold:
            break
    return less_p / n, s


def compute_pvalue_pickle(
    file_in,
    file_out,
    load_bar=False,
    S_min=2,
    S_max=5,
    mu_power=5,
    alpha_power=7,
    seed=None,
):
    data_out = dict()
    if seed is not None:
        np.random.seed(seed)  # Set seed for reproducibility in the main function
    with open(file_in, "rb") as f:
        data_in_pickle: dict = pickle.load(f)
    data_in = []
    for key, val in data_in_pickle.items():
        if np.any(np.isnan(val["r"])):
            data_out[key] = {
                "p_value": np.nan,
                "trials": np.nan,
            }
            continue
        data_in.append((key, val))

    # Calculate thresholds only once and pass them to inner functions
    thresholds = compute_thresholds(S_min, S_max, mu_power, alpha_power)

    with open(file_out, "wb") as f:
        pickle.dump(data_out, f)


if __name__ == "__main__":
    print("Running threshold with reproducible seed")
    t1 = time.time()
    compute_pvalue_pickle("input_file.pickle", "output_file.pickle", seed=42)
    t2 = time.time()
    print("Time elapsed:", t2 - t1)
