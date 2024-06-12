import numpy as np
import random
import time
from joblib import Parallel, delayed
from surface_code_utils import *

num_cores = 28
dist_list = [4,6,8]
p_bulk_qubit_list = np.linspace(1,10,8)*1e-2
Niter = 1000
Nrep = 28

params = {
    "p_bulk_qubit" : p_bulk_qubit_list,
    "p_bdy_qubit" : 0.1,
    "p_bulk_stab" : 0.01,
    "p_bdy_stab" : 0.1,
    "code_dist" : 4,
    "dx" : 4,
    "repetitions": 4,
    "Niter" : Niter,
    "loss_prob" : 0.1
}

loss_prob_list = [0.99]#np.linspace(0.1,0.8,5)

for loss_prob in loss_prob_list:
    params["loss_prob"] = loss_prob

    for i_d, d in enumerate(dist_list):
        dx = 2 * d
        dy = d
        repetitions = d
        code_details_2d = surface_code_star(dx, dy)
        code_details_3d = surface_code_3d_star(dx, dy,repetitions)
        code_details = code_details_2d, code_details_3d
        params["code_dist"] = d
        params["dx"] = dx
        params["repetitions"] = repetitions

        # logical_err_rate = joint_meas_2d_phase_diagram(code_details, params)
        # print(logical_err_rate)
        def runner(i_rep):
            tic = time.time()
            logical_err_rate = joint_meas_2d_phase_diagram(code_details, params)
            out_dir = "results/joint_meas_pd_star/"
            fname = f"d_{d}_star_pl_{loss_prob:.2f}_r_{i_rep}.npz"
            toc = time.time()
            print(f"{fname}, elapsed time {toc-tic:.2f} sec")
            np.savez(out_dir + fname, p_bulk_qubit_list, logical_err_rate)

        results = Parallel(n_jobs=num_cores)(delayed(runner)(i_rep) for i_rep in range(Nrep))

print("Job finished!")
