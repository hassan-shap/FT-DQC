import numpy as np
import random
import time
from joblib import Parallel, delayed
from surface_code_utils import *

num_cores = 28
dist_list = [4,6,8]
p_bulk_qubit_list = np.linspace(1,10,8)*1e-2
Niter = 1000
Nrep = 28*4

params = {
    "p_bulk_qubit" : p_bulk_qubit_list,
    "p_bdy_qubit" : 0.1,
    "p_bulk_stab" : 0.01,
    "p_bdy_stab" : 0.1,
    "code_dist" : 4,
    "dx" : 4,
    "repetitions": 4,
    "Niter" : Niter,
}

loss_prob_list = np.linspace(0.1,0.8,5)

for i_d, d in enumerate(dist_list):
    dx = 2 * d
    dy = d
    repetitions = d
    code_details = surface_code(dx, dy)
    params["code_dist"] = d
    params["dx"] = dx
    params["repetitions"] = repetitions

    # logical_err_rate = perf_sim_no_loss(code_details, params)
    # print(logical_err_rate)
    def runner(i_rep):
        tic = time.time()
        logical_err_rate = perf_sim_no_loss(code_details, params)
        out_dir = "results/joint_meas_pd/"
        fname = f"d_{d}_no_loss_r_{i_rep}.npz"
        toc = time.time()
        print(f"{fname}, elapsed time {toc-tic:.2f} sec")
        np.savez(out_dir + fname, p_bulk_qubit_list, logical_err_rate)

    results = Parallel(n_jobs=num_cores)(delayed(runner)(i_rep) for i_rep in range(Nrep))

print("Job finished!")
