import numpy as np
import time
from instance_maker import gen_instance_v3
from EM_simulate import simulate_Z
import json
import os

if __name__ == '__main__':
    # generate instance to test Z
    G_list = [2,4]
    I_list = [2,10]
    M = 50
    seed_list = [i for i in range(1,21)]
    J = 100
    step_size = 100
    samples = 10000
    for G in G_list:
        for I in I_list:
            for seed in seed_list:
                gen_instance_v3(G, I, M, J, name=f"instance_G{G}_I{I}_M{M}_J{J}_seed{seed}", 
                                terminar=False, seed = seed)

                # load instance with json instead of pickle
                with open(f"instances/instance_G{G}_I{I}_M{M}_J{J}_seed{seed}.json", 'r') as f:
                    data = json.load(f)
                X = np.array(data["n"])
                b = np.array(data["b"])
                p = np.array(data["p"])
                print(G,I)
                # simulate Z
                time_start = time.time()
                # create Z_instances folder if not there
                if not os.path.exists("Z_instances"):
                    os.makedirs("Z_instances")
                # Z = hit_and_run_matrix(X[0], p, b[0], G, I, samples = 1000, step_size=300, load_bar=True, verbose=True)
                Z = simulate_Z(X, p, b, samples = samples, step_size = step_size, verbose=True, load_bar=True, save=True, 
                               name=f"Z_instance_G{G}_I{I}_M{M}_J{J}_step{step_size}_S{samples}_seed{seed}_sequence", unique= False, parallel=True,
                               seed = seed)
                print(f"Time elapsed: {time.time() - time_start} seconds")

        # # save as pickle
    # with open(f"Z_instance_G{G}_I{I}_M{M}_J{J}_seed{0}.pickle", 'wb') as f:
    #     pickle.dump(Z, f)


    # # Z = [Z[s].to_list() for s in range(len(Z))]
    # # # save Z as json instead of pickle
    # # with open(f"Z_instance_G{G}_I{I}_M{M}_J{J}_seed{0}.json", 'w') as f:
    # #     json.dump(Z, f)

    # load instance with json instead of pickle

    # from instance_maker import gen_instance_v3

    # G,I,M,J = 10,10,5,100

    # gen_instance_v3(G,I,M,J, name=f"instance_G{G}_I{I}_M{M}_J{J}_seed{0}", terminar=False, seed = 0)

    # with open(f"instances/instance_G{G}_I{I}_M{M}_J{J}_seed{0}.json", 'r') as f:
    #     data = json.load(f)
    # X = np.array(data["n"])
    # b = np.array(data["b"])
    # p = np.array(data["p"])
    # import time
    # time_start = time.time()
    # Z = simulate_Z(X, p, b, samples = 100, step_size = 3000, verbose=True, load_bar=True, save=False, unique=True)
    # print(f"Time elapsed: {time.time() - time_start} seconds")


