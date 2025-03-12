import cvxpy as cp
import numpy as np
import pandas as pd

import pickle
import time

from pulp import LpMinimize, LpProblem, LpVariable, lpSum


from EM_mult import EM_mult
from helper_functions import get_partitions

"""Code for chosing optimal group aggregations"""

SOLVER_GUROBI = False


def std_group_combination(group_proportions, group_combination):
    # return np.std(np.sum(np.exp(group_proportions[:,group_combination]*50), axis=1))**2
    # return np.std(np.sum(-np.log(1-group_proportions[:,group_combination]), axis=1))**2
    return np.std(np.sum(group_proportions[:, group_combination], axis=1)) ** 2


# n is the number of age groups, groups should be ordered H (0,1,2,3,4,5,6,7) M (8,9,10,11,12,13,14,15)
def optimize_groups(group_matrix, age_groups=8, min_groups=1):

    ### GET SOME COMBINATORICS ###
    group_proportions = group_matrix / np.sum(group_matrix, axis=1)[:, None]
    a_H = [i for i in range(age_groups)]
    H = [
        a_H[i:j] for i in range(len(a_H)) for j in range(i + 1, len(a_H) + 1)
    ]  # TODAS LAS COMBINACIONES PARA EDADES HOMBRES
    a_M = [i for i in range(age_groups, 2 * age_groups)]
    M = [
        a_M[i:j] for i in range(len(a_M)) for j in range(i + 1, len(a_M) + 1)
    ]  # TODAS LAS COMBINACIONES PARA EDADES MUJER
    HM = [
        H[i] + M[i] for i in range(len(H))
    ]  # TODAS LAS COMBINACIONES PARA EDADES HOMBRES Y MUJERES
    group_combinations = (
        H + M + HM
    )  # TODAS LAS COMBINACIONES POSIBLES CON RESTRICCIONES
    std_grp_cmb = [
        std_group_combination(group_proportions, group_combination)
        for group_combination in group_combinations
    ]  # DESVIACION ESTANDAR PARA CADA COMBINACION

    #### OPTIMIZACION ####

    ## variables
    n = len(group_combinations)  # NUMERO DE COMBINACIONES
    v = std_grp_cmb  # std deviations
    x = cp.Variable(n, integer=True)  # 1 if group combination is selected, 0 otherwise
    y = cp.Variable(n)  # auxiliar variable to set the average in the objective function
    t = cp.Variable(1)  # 1/t = number of groups selected
    M = 1  # 2**(1000) # big M
    A = np.zeros((n, 16))  # A_ig = 1 if group g is in combination i, 0 otherwise
    for i in range(n):
        A[i, group_combinations[i]] = 1

    ## constraints ##
    constraints = []

    for g in range(16):
        # constraints.append(cp.sum(x@a[:,i]) == 1)
        constraints.append(
            cp.sum(y @ A[:, g]) == t
        )  # every group belongs to one combination

    constraints.append(cp.sum(y) == 1)  # to set average

    constraints.append(t >= 0)  # positive t

    constraints.append(t <= 1 / min_groups)  # ver si incluir más grupos

    constraints.append(y >= 0)  # positive y

    big_m = [
        y <= M * x,
        y - t <= M * (1 - x),
        y - t >= -M * (1 - x),
    ]  # big M constraints

    constraints += big_m

    ## objective ##
    objective = cp.Minimize(-cp.sum(v @ x))

    ## optimize ##
    prob = cp.Problem(objective, constraints)

    result = prob.solve()

    print("Average std deviation: ", -result)

    chosen_combs = np.where(x.value == 1)[0]

    final_combinations = [group_combinations[c] for c in chosen_combs]

    return final_combinations


# n is the number of age groups, groups should be ordered H (0,1,2,3,4,5,6,7) M (8,9,10,11,12,13,14,15)
# def optimize_groups_v2(group_matrix, age_groups=8, mu = 0.3, gamma = 0.2):
#     M_size, G_size = group_matrix.shape

#     ### GET SOME COMBINATORICS ###
#     group_proportions = group_matrix/np.sum(group_matrix, axis=1)[:,None]
#     a_H = [i for i in range(age_groups)]
#     H = [a_H[i:j] for i in range(len(a_H)) for j in range(i + 1, len(a_H) + 1)] # TODAS LAS COMBINACIONES PARA EDADES HOMBRES
#     a_M = [i for i in range(age_groups, 2*age_groups)]
#     M = [a_M[i:j] for i in range(len(a_M)) for j in range(i + 1, len(a_M) + 1)] # TODAS LAS COMBINACIONES PARA EDADES MUJER
#     HM = [H[i] + M[i] for i in range(len(H))] # TODAS LAS COMBINACIONES PARA EDADES HOMBRES Y MUJERES
#     group_combinations = H + M + HM # TODAS LAS COMBINACIONES POSIBLES CON RESTRICCIONES
#     std_grp_cmb = np.array([std_group_combination(group_proportions, group_combination) for group_combination in group_combinations]) # DESVIACION ESTANDAR PARA CADA COMBINACION

#     n = len(group_combinations) # NUMERO DE COMBINACIONES

#     A = np.zeros((n, 16)) # A_ig = 1 if group g is in combination i, 0 otherwise
#     for i in range(n):
#         A[i,group_combinations[i]] = 1
#     R = group_proportions.copy()
#     v = std_grp_cmb.copy() # std deviations
#     s = np.sqrt(v)
#     mean_R = np.mean(R, axis = 1) # revisar

#     M_ones = np.ones(shape=(M_size, 1))
#     cov = (group_proportions@A.T - M_ones @ mean_R.T @ A.T).T @ (R @ A.T - M_ones @ mean_R.T @ A.T)
#     corr = cov / (s[...,None] @ s[None,...])

#     #### OPTIMIZACION ####

#     ## variables


#     x = cp.Variable(n, integer = True) # 1 if group combination is selected, 0 otherwise
#     M = 1#2**(1000) # big M


#     ## constraints ##
#     constraints = []

#     for g in range(16):
#         #constraints.append(cp.sum(x@a[:,i]) == 1)
#         constraints.append(cp.sum(x@A[:,g]) == 1) # every group belongs to one combination

#     constraints.append(x <= std_grp_cmb / mu) # min variance

#     constraints.append(x >= 0) # positive x

#     for l in range(n):
#     	for h in range(n):
#     		if l != h:
#     			if corr[l,h] >= gamma:
#     				constraints.append(x[l] + x[h] <= 1)
#     ## objective ##

#     objective = cp.Minimize(-cp.sum(x))

#     ## optimize ##
#     prob = cp.Problem(objective, constraints)

#     result = prob.solve()

#     #print('Average std deviation: ', -result)

#     chosen_combs = np.where(x.value == 1)[0]

#     final_combinations = [group_combinations[c] for c in chosen_combs]

#     return final_combinations


#### GET NEW GROUPS ####
# Z_msgi
def new_group_Z(Z, combinations):
    new_Z = np.zeros((Z.shape[0], len(combinations), Z.shape[2]), dtype=int)
    for index, combination in enumerate(combinations):
        new_Z[:, index, :] = np.sum(Z[:, combination, :], axis=1)
    return np.unique(new_Z, axis=0)


def new_group_matrix(group_matrix, combinations, group_type="age_groups"):
    new_gm = np.zeros((group_matrix.shape[0], len(combinations)), dtype=int)
    new_names = []
    for index, combination in enumerate(combinations):
        new_gm[:, index] = new_group_votes(group_matrix, combination)
        if group_type == "age_groups":
            new_names.append(
                new_age_group_name(combination)
            )  # cambiar si hay otra estructura de nombres
        if group_type == "age_sex_groups":
            new_names.append(new_group_name(combination))

    return new_gm, new_names


def new_group_name(combination):
    min_age_list = np.array([10 * (i + 1) for i in range(8)] * 2)
    min_age_list[[0, 8]] = 18
    max_age_list = np.array([10 * (i + 2) - 1 for i in range(8)] * 2)
    # names_array = np.array(group_names)
    # new_group = names_array[combination]
    min_age = np.min(min_age_list[combination])
    max_age = np.max(max_age_list[combination])
    sex = "H" if np.min(combination) <= 7 else ""
    sex += "M" if np.max(combination) >= 8 else ""
    if max_age < 89:
        new_name = f"{min_age}-{max_age} {sex}"
    elif min_age == 18:
        new_name = f"{sex}"
    else:
        new_name = f"{min_age}+ {sex}"
    return new_name


def new_age_group_name(combination):
    min_age_list = np.array([10 * (i + 1) for i in range(8)])
    min_age_list[0] = 18
    max_age_list = np.array([10 * (i + 2) - 1 for i in range(8)])
    # names_array = np.array(group_names)
    # new_group = names_array[combination]
    min_age = np.min(min_age_list[combination])
    max_age = np.max(max_age_list[combination])
    if max_age < 89:
        new_name = f"{min_age}-{max_age}"
    else:
        new_name = f"{min_age}+"
    return new_name


def new_group_votes(group_matrix, combination):
    return np.sum(group_matrix[:, combination], axis=1)


def optimze_age_groups(
    group_matrix, age_groups=8, min_std=0.005, max_corr=0.25, add_corr=False
):
    M_size, G_size = group_matrix.shape

    group_proportions = group_matrix / np.sum(group_matrix, axis=1)[:, None]
    a_H = [i for i in range(age_groups)]
    group_combinations = [
        a_H[i:j] for i in range(len(a_H)) for j in range(i + 1, len(a_H) + 1)
    ]
    std_grp_cmb = np.array(
        [
            std_group_combination(group_proportions, group_combination)
            for group_combination in group_combinations
        ]
    )

    macro_group_proportions = np.zeros((M_size, len(group_combinations)))
    for i, g in enumerate(group_combinations):
        macro_group_proportions[:, i] = np.sum(group_proportions[:, g], axis=1)

    n = len(group_combinations)

    A = np.zeros((n, age_groups))
    for i in range(n):
        A[i, group_combinations[i]] = 1
    R = group_proportions.copy()
    v = std_grp_cmb.copy()
    v[age_groups - 1] = min_std
    s = np.sqrt(v)
    mean_R = np.mean(R, axis=1)

    M_ones = np.ones(shape=(M_size, 1))
    # print(group_proportions.shape, A.shape, M_ones.shape, mean_R.shape, std_grp_cmb.shape, s.shape)
    # print(group_proportions @ A.T)

    # cov = (group_proportions @ A.T - M_ones @ mean_R.T @ A.T).T @ (R @ A.T - M_ones @ mean_R.T @ A.T)
    cov = np.cov(macro_group_proportions.T)
    corr = cov / (s[..., None] @ s[None, ...])

    x = cp.Variable(n, integer=True)
    # if add_corr: y = cp.Variable((n,n), integer=True)
    if add_corr:
        y = [[cp.Variable() for j in range(i)] for i in range(n)]
    M = 1

    constraints = []
    for g in range(age_groups):
        constraints.append(cp.sum(x @ A[:, g]) == 1)

    constraints.append(x <= s / min_std)
    # constraints.append(x <= 1)
    constraints.append(x >= 0)

    for l in range(n):
        for h in range(l):
            if corr[l, h] >= max_corr:
                constraints.append(x[l] + x[h] <= 1)
            if True:
                constraints.append(x[l] + x[h] - y[l][h] <= 1)
                constraints.append(y[l][h] <= x[l])
                constraints.append(y[l][h] <= x[h])
            else:
                constraints.append(y[l][h] == 0)
        # if add_corr:
        #     for k in range(l, n):
        #         constraints.append(y[l*n + k] == 0)

    if add_corr:
        # objective = cp.Minimize(-cp.sum(x)*100 + cp.sum(cp.multiply(y,corr)))
        # objective = cp.Minimize(cp.sum(cp.multiply(y,corr)))
        objective = cp.Minimize(0)
        # objective = cp.Minimize(0)
    else:
        objective = cp.Minimize(-cp.sum(x) * 10 - cp.sum(x @ s))
        # objective = cp.Minimize(-cp.sum(x)*100)

    prob = cp.Problem(objective, constraints)

    result = prob.solve(verbose=True)
    # print(prob.status)
    # print(x.value)
    # print(prob.value)
    chosen_combs = np.where(x.value == 1)[0]

    final_combinations = [group_combinations[c] for c in chosen_combs]

    return final_combinations


# con un grupo forzar varianza infinito


def optimze_age_groups_v2(
    group_matrix,
    age_groups=8,
    min_std=0.005,
    max_corr=0.25,
    add_corr=False,
    add_std=False,
):
    M_size, G_size = group_matrix.shape

    group_proportions = group_matrix / np.sum(group_matrix, axis=1)[:, None]
    a_H = [i for i in range(age_groups)]
    group_combinations = [
        a_H[i:j] for i in range(len(a_H)) for j in range(i + 1, len(a_H) + 1)
    ]
    std_grp_cmb = np.array(
        [
            std_group_combination(group_proportions, group_combination)
            for group_combination in group_combinations
        ]
    )

    macro_group_proportions = np.zeros((M_size, len(group_combinations)))
    for i, g in enumerate(group_combinations):
        macro_group_proportions[:, i] = np.sum(group_proportions[:, g], axis=1)

    n = len(group_combinations)

    A = np.zeros((n, age_groups))
    for i in range(n):
        A[i, group_combinations[i]] = 1
    R = group_proportions.copy()
    v = std_grp_cmb.copy()
    v[age_groups - 1] = min_std
    s = np.sqrt(v)
    mean_R = np.mean(R, axis=1)

    M_ones = np.ones(shape=(M_size, 1))
    # print(group_proportions.shape, A.shape, M_ones.shape, mean_R.shape, std_grp_cmb.shape, s.shape)
    # print(group_proportions @ A.T)

    # cov = (group_proportions @ A.T - M_ones @ mean_R.T @ A.T).T @ (R @ A.T - M_ones @ mean_R.T @ A.T)
    # cov = np.cov(macro_group_proportions.T, bias=True)
    # corr = cov / (s[..., None] @ s[None, ...])
    corr = np.corrcoef(macro_group_proportions.T)
    corr[np.isnan(corr)] = 0
    # print("corr = ", corr)

    if SOLVER_GUROBI:
        x = cp.Variable(n, integer=True)
        # if add_corr: y = cp.Variable((n,n), integer=True)
        if add_corr:
            y = cp.Variable(n * n)
        M = 1

        constraints = []
        for g in range(age_groups):
            constraints.append(cp.sum(x @ A[:, g]) == 1)

        constraints.append(x <= s / min_std)
        # constraints.append(x <= 1)
        constraints.append(x >= 0)

        for l in range(n):
            for h in range(l):
                if corr[l, h] >= max_corr:
                    constraints.append(x[l] + x[h] <= 1)
                if add_corr:
                    constraints.append(x[l] + x[h] - y[l * n + h] <= 1)
                    constraints.append(y[l * n + h] <= x[l])
                    constraints.append(y[l * n + h] <= x[h])
                else:
                    constraints.append(y[l * n + h] == 0)
            # if add_corr:
            #     for k in range(l, n):
            #         constraints.append(y[l*n + k] == 0)

        if add_corr:
            # objective = cp.Minimize(-cp.sum(x)*100 + cp.sum(cp.multiply(y,corr)))
            # objective = cp.Minimize(cp.sum(cp.multiply(y,corr)))
            objective = cp.Minimize(0)
            # objective = cp.Minimize(0)
        else:
            objective = cp.Minimize(-cp.sum(x) * 10 - cp.sum(x @ s))
            # objective = cp.Minimize(-cp.sum(x)*100)

        prob = cp.Problem(objective, constraints)

        result = prob.solve(verbose=True)
        # print(prob.status)
        # print(x.value)
        # print(prob.value)
        chosen_combs = np.where(x.value == 1)[0]

    ######### NOT GUROBI #### THIS WORKS
    else:
        model = LpProblem("LP", LpMinimize)
        x = [LpVariable(f"x_{i}", cat="Binary") for i in range(n)]
        if add_corr:
            y = [
                [LpVariable(f"y_{i}{j}", lowBound=0, upBound=1) for j in range(i)]
                for i in range(n)
            ]
        for g in range(age_groups):
            model += lpSum(x @ A[:, g]) == 1
        for i in range(n):
            model += x[i] <= s[i] / min_std
        for l in range(n):
            for h in range(l):
                # for h in range(n):
                if corr[l, h] >= max_corr:
                    model += x[l] + x[h] <= 1
                if add_corr:
                    model += x[l] + x[h] - y[l][h] <= 1
                    model += y[l][h] <= x[l]
                    model += y[l][h] <= x[h]
                # else:
                #     model += y[l][h] == 0
            # if add_corr:
            #     for k in range(l, n):
            #         model += y[l][k] <= 0
        if add_corr:
            objective = -10 * lpSum(x) + lpSum(
                [y[i][j] * corr[i, j] for i in range(n) for j in range(i)]
            )
        elif add_std:
            objective = -lpSum(x) * 10 - lpSum(x @ s)
        else:
            objective = -lpSum(x)

        model += objective
        # pulp.GUROBI(msg=0).solve(model)
        model.solve()
        xx = [x[i].value() for i in range(n)]
        # yy = [[y[i][j].value() for j in range(i)] for i in range(n)]
        # print("xx = ", xx)
        # print("yy = ", yy)
        chosen_combs = np.where(np.array(xx) == 1)[0]
        # print("chosen_combs = ", chosen_combs)

    final_combinations = [group_combinations[c] for c in chosen_combs]

    return final_combinations


# study of the variance of the groups
def bootstrap(x, b, S=20, seed=None):
    """returns the probability matrix for S bootstraps"""
    # choose random index for ballot-boxes
    np.random.seed(seed)
    M = len(x)
    G = b.shape[1]
    I = x.shape[1]
    p_bootstrap = np.zeros((S, G, I))
    for s in range(S):
        idx = np.random.choice(M, size=(M), replace=True)
        x_ = x[idx, :].copy()
        b_ = b[idx, :].copy()
        results = EM_mult(x_, b_, max_iterations=1000, verbose=False, load_bar=False)
        if np.isnan(results[0]).any():
            print("x", x_)
            print("b", b_)
        p_bootstrap[s] = results[0]
    return p_bootstrap


def group_election_routine(
    x, b, std_thres=0.02, age_groups=8, threshold="g", S=20, log_func=False
):
    """returns the biggest group combinations that have a std deviation below std_thres
    returns the prbability std deviation of the biggest group combination
    """
    # print('trying', 8)
    p_bootstrap = bootstrap(x, b, S=S)
    if threshold == "gc":
        p_std = np.std(p_bootstrap, axis=0)
    if threshold == "g":
        p_std = np.mean(np.std(p_bootstrap, axis=0), axis=1)
    if threshold == "p":
        p_std = np.mean(np.std(p_bootstrap, axis=0))
    # print(p_std)
    # print(p_std)
    if np.all(p_std < std_thres):
        return [[i] for i in range(age_groups)], p_std
    else:
        group_matrix = b.copy()
        M_size, G_size = group_matrix.shape

        group_proportions = group_matrix / np.sum(group_matrix, axis=1)[:, None]
        a_H = [i for i in range(age_groups)]
        group_combinations = [
            a_H[i:j] for i in range(len(a_H)) for j in range(i + 1, len(a_H) + 1)
        ]
        std_grp_cmb = np.array(
            [
                std_group_combination(group_proportions, group_combination)
                for group_combination in group_combinations
            ]
        )

        macro_group_proportions = np.zeros((M_size, len(group_combinations)))
        for i, g in enumerate(group_combinations):
            macro_group_proportions[:, i] = np.sum(group_proportions[:, g], axis=1)

        n = len(group_combinations)

        A = np.zeros((n, age_groups))
        for i in range(n):
            A[i, group_combinations[i]] = 1
        R = group_proportions.copy()
        v = std_grp_cmb.copy()
        # v[age_groups - 1] = 10
        s = np.sqrt(v)
        if log_func:
            s = np.log(s)
        for g_size in range(7, 1, -1):
            # print('trying', g_size)
            XG = cp.Variable(n, boolean=True)
            constraints = []
            for g in range(age_groups):
                constraints.append(cp.sum(XG @ A[:, g]) == 1)
            constraints.append(cp.sum(XG) == g_size)
            objective = cp.Minimize(-cp.sum(XG @ s))
            prob = cp.Problem(objective, constraints)
            result = prob.solve(verbose=False)
            chosen_combs = np.where(XG.value == 1)[0]
            final_combinations = [group_combinations[c] for c in chosen_combs]

            b_new, _ = new_group_matrix(b, final_combinations)
            results = {}
            p_bootstrap = bootstrap(x, b_new, S=20, seed=123)
            if threshold == "gc":
                p_std = np.std(p_bootstrap, axis=0)
            if threshold == "g":
                p_std = np.mean(np.std(p_bootstrap, axis=0), axis=1)
            if threshold == "p":
                p_std = np.mean(np.std(p_bootstrap, axis=0))
            # print(p_std)
            if np.all(p_std < std_thres):
                return final_combinations, p_std
        final_combinations = [[i for i in range(age_groups)]]
        b_new, _ = new_group_matrix(b, final_combinations)
        p_bootstrap = bootstrap(x, b_new, S=S, seed=123)
        if threshold == "gc":
            p_std = np.std(p_bootstrap, axis=0)
        if threshold == "g":
            p_std = np.mean(np.std(p_bootstrap, axis=0), axis=1)
        if threshold == "p":
            p_std = np.mean(np.std(p_bootstrap, axis=0))
        return final_combinations, p_std


# old version
def evaluate_all_combs_A_old(x, b, A, S=20, std_thres=0.05):
    # get all possible combinations of size A [0,1,2,3,4,5,6,7] (place A-1 separations)
    G = b.shape[1]
    group_combs_A = get_partitions(G, A)
    combination_list = []
    ll_list = []
    under_threshold = []
    obj_func_list = []
    log_obj_func_list = []
    p_std_max_list = []
    for combination in group_combs_A:
        print(combination)
        b_new, _ = new_group_matrix(b, combination)
        b_prop = b_new / np.sum(b_new, axis=1)[..., None]  # get proportions (normalize)
        b_prop_std = np.std(b_prop, axis=0)
        obj_func = np.sum(b_prop_std)
        log_obj_func = np.sum(np.log(b_prop_std))
        p_bootstrap = bootstrap(x, b_new, S=S)
        p_std = np.std(p_bootstrap, axis=0)

        # compute EM_algorithm
        results = {}
        em = EM_mult(
            x,
            b_new,
            max_iterations=10000,
            dict_results=results,
            verbose=False,
            load_bar=False,
        )
        ll = results["ll"]
        if np.all(p_std < std_thres):
            under_threshold.append(True)
        else:
            under_threshold.append(False)
        print(combination, ll, under_threshold[-1])
        combination_list.append(combination)
        ll_list.append(ll)
        obj_func_list.append(obj_func)
        log_obj_func_list.append(log_obj_func)
        p_std_max_list.append(np.max(p_std))
    return (
        combination_list,
        ll_list,
        under_threshold,
        obj_func_list,
        log_obj_func_list,
        p_std_max_list,
    )


# new version: store all results in a dictionary
def evaluate_all_combs_A(x, b, A, S=20, std_thres=0.05, results_dict={}, circ_name=""):
    G = b.shape[1]
    group_combs_A = get_partitions(G, A)
    # results_dict = {}  # Store all results in a dictionary

    for combination in group_combs_A:
        print(combination)
        b_new, _ = new_group_matrix(b, combination)
        b_prop = b_new / np.sum(b_new, axis=1)[..., None]  # get proportions (normalize)
        b_prop_std = np.std(b_prop, axis=0)
        b_std = np.std(b_new, axis=0)
        obj_func_prop = np.sum(b_prop_std)
        log_obj_func_prop = np.sum(np.log(b_prop_std))
        obj_func = np.sum(b_std)
        log_obj_func = np.sum(np.log(b_std))
        p_bootstrap = bootstrap(x, b_new, S=S)
        print("The bootstrap is")
        p_std = np.std(p_bootstrap, axis=0)
        print(p_std)
        print("The mean is:", np.mean(p_std))

        # Compute EM algorithm
        em_results = {}
        EM_mult(
            x,
            b_new,
            max_iterations=10000,
            dict_results=em_results,
            verbose=False,
            load_bar=False,
        )
        ll = em_results["ll"]

        combination_name = ""
        for gr in combination:
            combination_name += new_age_group_name(gr) + ","
        # remove last character
        combination_name = combination_name[:-1]

        # Store results in the dictionary using the combination as key
        results_dict[(combination_name, circ_name)] = {
            "circ": circ_name,
            "M": x.shape[0],
            "A": A,
            "comb name": combination_name,
            "log_likelihood": ll,
            "under_threshold": np.all(p_std < std_thres),
            "obj_func": obj_func,
            "log_obj_func": log_obj_func,
            "obj_func_prop": obj_func_prop,
            "log_obj_func_prop": log_obj_func_prop,
            "p_std_max": np.max(p_std),
            "combination": combination,
        }

        # print(combination, ll, results_dict[tuple(combination)]['under_threshold'])

    # return results_dict


if __name__ == "__main__":
    """
    TODO:
    - correr distintas obj con más instancias (Pablo)
    """
    n_combs = 50
    electors = pd.read_csv(
        "2021_11_Presidencial/output/2021_11_Presidencial_ELECTORES.csv"
    )
    votes = pd.read_csv("2021_11_Presidencial/output/2021_11_Presidencial_VOTOS.csv")
    # READ PICKLE
    with open("2021_11_Presidencial/output/GRUPOS.pickle", "rb") as handle:
        grupos = pickle.load(handle)
    with open("2021_11_Presidencial/output/CANDIDATOS.pickle", "rb") as handle:
        candidatos = pickle.load(handle)

    llave_mesa = ["CIRCUNSCRIPCION ELECTORAL", "MESA"]
    df = votes.merge(electors[llave_mesa + grupos], on=llave_mesa, how="inner")

    circs = electors["CIRCUNSCRIPCION ELECTORAL"].unique()

    # select 10 at random
    # np.random.seed(123)
    np.random.seed(n_combs)
    # circs = np.random.choice(circs, size = n_combs, replace = False)

    # MIP_df = pd.read_csv('MIP_df_v4.csv')
    # circs = ['ORILLA DEL MAULE', 'VITACURA', 'MAIPU']
    # circs = ['CALAMA']
    # circs += list(MIP_df['Circunscripcion'].unique())

    circs = ["NIEBLA"]
    wrong_circs = []
    results = []

    results_dict = {}
    for c in circs:
        print("-" * 40)
        print(c)
        df_c = df.loc[df["CIRCUNSCRIPCION ELECTORAL"] == c]
        x = df_c[candidatos].values
        b = df_c[grupos].values
        # start = time.time()
        # optimal_comb, p_std = group_election_routine(x, b, std_thres = 0.05, age_groups = 8, S=20, threshold = 'gc')
        # end = time.time()
        # results.append([c, x.shape[0] ,len(optimal_comb), optimal_comb, p_std, end-start])
        # print('Selected combination', optimal_comb)
        # # print(c, x.shape[0], combinations, end-start)
        # # print time it took
        # print(f'time routine {end-start}')

        # A = len(optimal_comb)
        A = 1
        G = 8
        start = time.time()
        # try all combinations
        for a in range(G, A - 1, -1):
            print(f"evaluating all combs size {a}")
            np.random.seed(1)
            evaluate_all_combs_A(
                x, b, a, S=400, std_thres=0.0001, results_dict=results_dict, circ_name=c
            )
            # combination_list, ll_list, under_threshold, obj_function_list, log_obj_function_list, p_std_list = evaluate_all_combs_A(x, b, a, S=20, std_thres = 0.05)
            # optimal = False
            # for ix in range(len(combination_list)):
            #     # optimal = False
            #     # if combination_list[ix] == optimal_comb:
            #     #     optimal = True
            #     # combination name
            #     combination_name = ''
            #     for gr in combination_list[ix]:
            #         combination_name += new_age_group_name(gr) + ','
            #     # remove last character
            #     combination_name = combination_name[:-1]
            #     # results.append([c, combination_list[ix], combination_name, ll_list[ix], under_threshold[ix], optimal, a])
            #     results.append([c, a, x.shape[0], combination_list[ix], combination_name, ll_list[ix],
            #                     under_threshold[ix], obj_function_list[ix], log_obj_function_list[ix],
            #                     p_std_list[ix]])

            # if a == A:
            # look for optimal

            # try:
        print(f"time all combs {time.time()-start}")
        # except:
        #     wrong_circs.append((c, electors_c.shape[0]))
    print(len(circs))

    # save results
    # results_df = pd.DataFrame(results, columns = ['circ', 'A', 'M', 'comb', 'comb_name', 'll', 'under_threshold', 'obj_func', 'log_obj_func', 'p_std'])
    results_df = pd.DataFrame.from_dict(results_dict, orient="index")
    # results_df.to_csv(f'2021_11_Presidencial/combinations/all_combinations_{n_combs}.csv', index = False)

    # results_df = pd.DataFrame(results, columns = ['circ', 'combination', 'combination_name', 'll', 'under_threshold', 'MIP_optimal', 'A'])
    # results_df.to_csv('2021_11_Presidencial/combinations/combinations.csv', index = False)

    # results_df = pd.DataFrame(results, columns = ['circ', 'M', 'G', 'combinations', 'p_std', 'time'])
    # results_df.to_csv('group_routine_results_v2.csv', index = False)

    # test MIP
    # electors = pd.read_csv('2021_11_Presidencial/output/2021_11_Presidencial_ELECTORES.csv')
    # votes = pd.read_csv('2021_11_Presidencial/output/2021_11_Presidencial_VOTOS.csv')
    # # READ PICKLE
    # with open('2021_11_Presidencial/output/GRUPOS.pickle', 'rb') as handle:
    #     grupos = pickle.load(handle)

    # with open('2021_11_Presidencial/output/CANDIDATOS.pickle', 'rb') as handle:
    #     candidatos = pickle.load(handle)
    # llave_mesa = ['CIRCUNSCRIPCION ELECTORAL', 'MESA']
    # df = votes.merge(electors[llave_mesa + grupos], on = llave_mesa, how = 'inner')

    # np.random.seed(123)
    # # count how many rows each "CIRCUNSCRIPCION ELECTORAL" has
    # num_mesas = votes.groupby('CIRCUNSCRIPCION ELECTORAL').count()['MESA']
    # rango_mesas = [[3,10],[10,50],[50,100],[100,200],[200,300],[400,600]]
    # S_mesas = 5
    # circs = []
    # for m in rango_mesas:
    #     circ_m = num_mesas[(num_mesas >= m[0] ) & (num_mesas < m[1])]
    #     # choose random index
    #     idx = np.random.choice(circ_m.index, size = S_mesas, replace = False)
    #     circs += list(idx)

    # min_stds = [0.01, 0.05, 0.1, 0.2, 0.3]
    # # min_stds = [0.1]
    # # max_corrs = [0, 0.2, 0.4, 0.6, 0.8, 1]
    # max_corrs = [0.6, 0.8, 1]

    # tie = ['corr', 'std', 'none']

    # MIP_df = []

    # for circ in circs:
    #     print(circ)
    #     opt_ll = -np.inf
    #     x = df.loc[df['CIRCUNSCRIPCION ELECTORAL'] == circ, candidatos].values

    #     # # JUNTAR VOTOS BLANCOS Y NULOS
    #     # x[:,6] = x[:,6] + x[:,7]
    #     # x = x[:,[0,1,2,3,4,5,6,8]]

    #     b = df.loc[df['CIRCUNSCRIPCION ELECTORAL'] == circ, grupos].values
    #     M = len(x)
    #     for min_std in min_stds:
    #         for max_corr in max_corrs:
    #             for i, t in enumerate(tie):
    #                 print(min_std, max_corr, t)
    #                 if i == 0:
    #                     group_combinations = optimze_age_groups_v2(b, age_groups = 8, min_std = min_std, max_corr = max_corr, add_corr = True)
    #                 if i == 1:
    #                     group_combinations = optimze_age_groups_v2(b, age_groups = 8, min_std = min_std, max_corr = max_corr, add_std = True)
    #                 if i == 2:
    #                     group_combinations = optimze_age_groups_v2(b, age_groups = 8, min_std = min_std, max_corr = max_corr)
    #                 b_new, _ = new_group_matrix(b, group_combinations)
    #                 proportion = b_new/np.sum(b_new, axis = 1)[...,None]
    #                 # print(np.corrcoef(proportion.T))
    #                 results = {}
    #                 em = EM_mult(x, b_new, max_iterations = 10000, dict_results = results, verbose = False, load_bar =False)
    #                 ll = results['Q'] - results['E_log_q']
    #                 p_bootstrap = bootstrap(x, b_new, S = 100, seed = 123)
    #                 p_std = np.mean(np.std(p_bootstrap, axis = 0))
    #                 MIP_df.append([circ, M, min_std, max_corr, group_combinations, ll, p_std, t])

    # MIP_df = pd.DataFrame(MIP_df, columns = ['Circunscripcion', 'M', 'min_std', 'max_corr', 'group_combinations', 'll', 'p_std', 'tie'])
    # MIP_df.to_csv('MIP_df_v4.csv')

    # print(wrong_circs)
