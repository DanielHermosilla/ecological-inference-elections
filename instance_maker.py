import numpy as np
import random
import pickle
import json
import ast
import os

"""Code for generating election instances"""


def gen_instance(M, J, p_gi, lambda_, name=None, terminar=False, seed=None):
    # screen_print = input("Deseas imprimir el resultado en pantalla? (y/n)\n")
    # while screen_print not in ("y", "n"):
    #     screen_print = input("Por favor escribe \"y\" para si o \"n\" para no\n")
    # if screen_print == "y":
    #     screen_print = True
    # else:
    #     screen_print = False

    # save_file = input("Deseas guardar el resultado en un archivo? (y/n)\n")
    # while save_file not in ("y", "n"):
    #     save_file = input("Por favor escribe \"y\" para si o \"n\" para no\n")
    # if save_file == "y":
    #     save_file = True
    # else:
    #     save_file = False

    G = len(p_gi)
    I = len(p_gi[0])
    N = M * J
    np.random.seed(seed)
    G_fractions = np.random.dirichlet([10.0] * G, size=None)
    G_population = np.repeat(0, N)
    k = 1
    for i in range(N):
        G_population[i] = k
        if i / N > sum(G_fractions[0:k]):
            k += 1
    indices_aleatorizar = random.sample(range(N), round(lambda_ * N))
    valores_aleatorizar = G_population[indices_aleatorizar]

    indices_aleatorizar.sort()
    for k in range(len(valores_aleatorizar)):
        G_population[indices_aleatorizar[k]] = valores_aleatorizar[k]
    G_population_2 = G_population.reshape(M, J)
    b_mg = [[sum(G_population_2[m] == g + 1) for g in range(0, G)] for m in range(0, M)]
    r_mgi = [
        [np.random.multinomial(b_mg[m][g], p_gi[g]).tolist() for g in range(0, G)]
        for m in range(0, M)
    ]

    n_mi = n_mi_from_r_mgi(r_mgi)

    # if(save_file):
    #    print("IMPORTANTE: Si escribes un nombre que ya existe el archivo se sobreescribira!")
    data = {}
    data["M"] = M
    data["G"] = G
    data["I"] = I
    data["J"] = J
    data["p"] = np.array(p_gi)
    data["r"] = np.array(r_mgi)
    data["n"] = np.array(n_mi)
    data["lambda_"] = lambda_
    aux_b_mg = []
    for i in b_mg:
        aux_aux = []
        for j in i:
            aux_aux.append(int(j))
        aux_b_mg.append(aux_aux)
    data["b"] = aux_b_mg
    if name is None:
        name = input("Que nombre deseas para tu archivo?\n")
    # save as pickle
    with open(f"instances/{name}.pickle", "wb") as f:
        pickle.dump(data, f)

    # with open(f"forced_instances/{file_name}.json", 'w') as f:
    #     json.dump(data, f)
    if terminar:
        exit()


def gen_instance_v2(
    G, I, M, J, p_gi=None, lambda_=0.5, name=None, terminar=False, seed=None
):
    np.random.seed(seed)
    if p_gi is None:
        p_gi = np.random.dirichlet([1] * I, size=G)

    N = M * J

    # G_fractions = np.random.dirichlet([10.0]*G, size=None) # GRUPOS ALERATORIOS
    G_fractions = [1 / G] * G  # TODOS LOS GRUPOS TIENEN LA MISMA CANTIDAD DE VOTANTES
    G_population = np.repeat(0, N)
    k = 1
    for i in range(N):
        G_population[i] = k
        if i / N > sum(G_fractions[0:k]):
            k += 1
    # choose random index with numpy
    indices_aleatorizar = np.random.choice(N, round(lambda_ * N), replace=False)
    # indices_aleatorizar = random.sample(range(N), round(lambda_ * N))
    valores_aleatorizar = G_population[indices_aleatorizar]

    indices_aleatorizar.sort()
    for k in range(len(valores_aleatorizar)):
        G_population[indices_aleatorizar[k]] = valores_aleatorizar[k]
    G_population_2 = G_population.reshape(M, J)
    b_mg = [
        [int(sum(G_population_2[m] == g + 1)) for g in range(0, G)] for m in range(0, M)
    ]
    r_mgi = [
        [
            np.random.multinomial(b_mg[m][g], p_gi[g]).astype(int).tolist()
            for g in range(0, G)
        ]
        for m in range(0, M)
    ]

    n_mi = n_mi_from_r_mgi(r_mgi)
    # if(save_file):
    #    print("IMPORTANTE: Si escribes un nombre que ya existe el archivo se sobreescribira!")
    data = {}
    data["M"] = M
    data["G"] = G
    data["I"] = I
    data["J"] = J
    data["p"] = p_gi.tolist()
    data["r"] = r_mgi
    data["n"] = n_mi
    data["lambda_"] = lambda_
    data["b"] = b_mg
    if name is None:
        name = input("Que nombre deseas para tu archivo?\n")
    # save as pickle
    # with open(f"instances/{name}.pickle", 'wb') as f:
    #     pickle.dump(data, f)
    with open(f"instances/{name}.json", "w") as f:
        json.dump(data, f)
    if terminar:
        exit()


def gen_instance_v3(
    G, I, M, J, p_gi=None, lambda_=0.5, name=None, terminar=False, seed=None
):
    np.random.seed(seed)
    if p_gi is None:
        p_gi = np.random.dirichlet([1] * I, size=G)

    N = M * J

    # G_fractions = np.random.dirichlet([10.0]*G, size=None) # GRUPOS ALERATORIOS
    G_fractions = [1 / G] * G  # TODOS LOS GRUPOS TIENEN LA MISMA CANTIDAD DE VOTANTES
    G_population = np.repeat(0, N)
    votes_population = np.repeat(0, N)
    k = 1
    for i in range(N):
        G_population[i] = k
        votes_population[i] = np.random.choice(I, p=p_gi[k - 1])
        if i / N > sum(G_fractions[0:k]):
            k += 1
    # choose random index with numpy
    indices_aleatorizar = np.random.choice(N, round(lambda_ * N), replace=False)
    # indices_aleatorizar = random.sample(range(N), round(lambda_ * N))
    valores_aleatorizar = G_population[indices_aleatorizar]
    votos_aleatorizar = votes_population[indices_aleatorizar]
    indices_aleatorizar.sort()
    for k in range(len(valores_aleatorizar)):
        G_population[indices_aleatorizar[k]] = valores_aleatorizar[k]
        votes_population[indices_aleatorizar[k]] = votos_aleatorizar[k]

    G_population_2 = G_population.reshape(M, J)
    votes_population_2 = votes_population.reshape(M, J)
    b_mg = [
        [int(sum(G_population_2[m] == g + 1)) for g in range(0, G)] for m in range(0, M)
    ]
    r_mgi = [
        [
            [
                int(sum(votes_population_2[m][G_population_2[m] == g + 1] == i))
                for i in range(I)
            ]
            for g in range(0, G)
        ]
        for m in range(0, M)
    ]
    # r_mgi = [[np.random.multinomial(b_mg[m][g], p_gi[g]).astype(int).tolist() for g in range(0, G)] for m in range(0, M)]
    n_mi = n_mi_from_r_mgi(r_mgi)

    # if(save_file):
    #    print("IMPORTANTE: Si escribes un nombre que ya existe el archivo se sobreescribira!")
    data = {}
    data["M"] = M
    data["G"] = G
    data["I"] = I
    data["J"] = J
    data["p"] = p_gi.tolist()
    data["r"] = r_mgi
    data["n"] = n_mi
    data["lambda_"] = lambda_
    data["b"] = b_mg
    if name is None:
        name = input("Que nombre deseas para tu archivo?\n")
    # save as pickle
    # with open(f"instances/{name}.pickle", 'wb') as f:
    #     pickle.dump(data, f)
    with open(f"instances/{name}.json", "w") as f:
        json.dump(data, f)
    if terminar:
        exit()


# def gen_instance_dirichlet(M, J, G, C, p_gi, lambda_, name = None, terminar = False, seed = None):


def load_instance():
    nombre_instancia = input("Ingrese nombre instancia a cargar.\n")
    # with open(os.getcwd()+f"/instances/{nombre_instancia}.pickle", 'rb') as file:
    #     data = pickle.load(file)
    # read json file
    with open(os.getcwd() + f"/instances/{nombre_instancia}.json", "r") as file:
        data = json.load(file)
    return data
    # M = data["M"]
    # J = data["J"]
    # p_gi = data["p_gi"]
    # lambda_ = data["lambda_"]
    # r_mgi = data["r_mgi"]
    # n_mi = data["n_mi"]
    # b_mg = data["b_mg"]
    # return M, J, p_gi, lambda_, r_mgi, n_mi, b_mg


def n_mi_from_r_mgi(r_mgi):
    M = len(r_mgi)
    G = len(r_mgi[0])
    I = len(r_mgi[0][0])
    n_mi = [
        [sum([r_mgi[m][g][i] for g in range(0, G)]) for i in range(0, I)]
        for m in range(0, M)
    ]
    return n_mi


def print_r_mgi_b_gi(r_mgi, b_mg):

    M = len(r_mgi)
    G = len(r_mgi[0])
    I = len(r_mgi[0][0])
    print("\n\ndef force_instance(printear_instancia = False):")
    print("\tM = " + str(M))
    print("\tG = " + str(G))
    print("\tI = " + str(I))
    print("\tr_mgi = [[[0 for i in range(I)] for g in range(G)] for m in range(M)]")
    print("\tb_mg = [[0 for g in range(G)] for m in range(M)]\n")
    for m in range(0, M):
        for g in range(0, G):
            for i in range(0, I):
                print(
                    "\tr_mgi["
                    + str(m)
                    + "]["
                    + str(g)
                    + "]["
                    + str(i)
                    + "] = "
                    + str(r_mgi[m][g][i])
                )
    print("")
    for m in range(0, M):
        for g in range(0, G):
            print("\tb_mg[" + str(m) + "][" + str(g) + "] = " + str(b_mg[m][g]))
    print("\n\tn_mi = n_mi_from_r_mgi(r_mgi)")
    print("\tif printear_instancia:")
    print("\t\tprint_instance(r_mgi, b_mg, n_mi)")
    print("\treturn r_mgi, n_mi, b_mg\n")


if __name__ == "__main__":
    print("Bienvenido a Instance Maker!")
    I = int(input("Por favor coloca el valor de I (candidatos) \n I:"))
    G = int(input("Por favor coloca el valor de G (grupos) \n G:"))
    M = int(input("Por favor coloca el valor de M (mesas) \n M:"))
    J = int(input("Por favor coloca el valor de J (personas) \n J:"))
    # I = 2
    # G = 2
    # M = 10
    # J = 100
    p_gi = np.array([[0.4, 0.6], [0.7, 0.3]])

    p_gi = np.zeros((G, I))
    for g in range(G):
        for i in range(I):
            prob = input(
                "Por favor coloca el valor de p_gi (probabilidades voto grupo-candidato) \n p_gi:"
            )
            p_gi[g][i] = prob
    # I = int(input("Por favor coloca el valor de I (candidatos) \n I:"))
    # G = int(input("Por favor coloca el valor de G (grupos) \n G:"))
    # M = int(input("Por favor coloca el valor de M (mesas) \n M:"))
    # J = int(input("Por favor coloca el valor de J (personas) \n J:"))
    I = 2
    G = 2
    M = 10
    J = 100
    p_gi = np.array([[0.4, 0.6], [0.7, 0.3]])

    # p_gi = input("Por favor coloca el valor de p_gi (probabilidades voto grupo-candidato) \n p_gi:")
    # lambda_ = float(input("Por favor coloca el valor de lambda_ (distribucion votantes) \n lambda_:"))
    # p_gi = ast.literal_eval(p_gi)
    # p_gi = [[n for n in p] for p in p_gi]
    # gen_instance_v2(G, I, M, J, p_gi=p_gi, name="A1", terminar=False, seed=None)
    gen_instance_v3(G, I, M, J, p_gi=p_gi, name="A2", terminar=False, seed=None)
