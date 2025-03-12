import numpy as np
import time
from tqdm import tqdm
from EM_full import EM_full
from EM_mult import EM_mult

def compute_q_binom(n,p,b):
    M, G = b.shape
    I = p.shape[1]
    q = np.zeros(shape=(M,G,I))
    for m in range(M):
        q[m] = compute_qm_binom(n[m], p, b[m], I, G)
    #q[np.isnan(q)] = 0
    return q

def compute_qm_binom(n_m, p, b_m, I, G):
    q_m = np.zeros(shape=(G,I))
    J = np.sum(n_m)
    for g in range(G):
        for i in range(I):
            q_m[g,i] = np.math.factorial(b_m[g])*(p[g,i]**(n_m[i] - 1) * (1-p[g,i])**(b_m[g]-n_m[i]))/np.math.factorial(n_m[i]-1)
    q_m = q_m/np.sum(q_m, axis=1)[:,None]
    return q_m

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
    num = np.sum(np.multiply(q,b[...,None]),axis=0)
    dem = np.sum(b,axis=0)[...,None]
    return num / dem

# (g,i) compute estimate of p using EM algorithm with parameters X and b 
def EM_binom(X, b, convergence_value = 0.0001, max_iterations = 100, 
                 p_method = 'group_proportional', load_bar = True, verbose = True):
    M_size, I_size = X.shape
    G_size = b.shape[1] 
    J_mean = np.round(np.mean(np.sum(X, axis = 1)),1)
    if verbose: print('M =',M_size,' G =',G_size,' I =',I_size,' J =',J_mean,' delta =',convergence_value)
    if p_method == 'uniform':
        p_est = np.full((G_size, I_size), 1/I_size)
    if p_method == 'proportional':
        p_est = np.array([np.sum(X, axis = 0) / np.sum(X) for g in range(G_size)])
    if p_method == 'group_proportional':
        p_mesas = X/np.sum(X,axis=1)[...,None]
        p_est = np.zeros((G_size, I_size))
        for g in range(G_size):
            for i in range(I_size):
                p_est[g,i] = np.sum(p_mesas[:,i]*b[:,g])/np.sum(b[:,g])
        
    if verbose:
        print('-'*100)
        print('EM-algorithm')
        print('-'*100)

    start = time.time()
    for i in tqdm(range(1,max_iterations+1), disable = not load_bar):
        q = compute_q_binom(X, p_est, b)
        p_new = compute_p(q,b)
        if verbose:
            print('-'*50)
            print(p_est)
            print('Δ: ', np.max(np.abs(p_new-p_est)))
            print('-'*50)
        if (np.abs(p_new-p_est) < convergence_value).all():
            end = time.time()
            if verbose: print(f'Convergence took {i} iterations and {end-start} seconds.')
            return p_new, i, end-start
        p_est = p_new.copy()
    end = time.time()
    if verbose: print(f'Did not reach convergence after {i} iterations and {end-start} seconds. \n')
    return p_est, i, end-start


if __name__ == '__main__':
    # load json
    import json
    with open('instances\instance_G2_I2_M50_J50.json') as json_file:
        data = json.load(json_file)       
    X = np.array(data['n'])
    b = np.array(data['b'])
    p = np.array(data['p'])
    print(p)
    print(EM_binom(X, b, p_method = 'group_proportional', verbose = False))
    print(EM_full(X, b, p_method = 'group_proportional', verbose = False))
    print(EM_mult(X, b, p_method = 'group_proportional', verbose = False))