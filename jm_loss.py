import numpy as np
import random
import time
from joblib import Parallel, delayed
from surface_code_utils import *

num_cores = 28
dist_list = [4,6,8]
bandwidth = 4
gen_coh_ratio_list = np.logspace(-3,-1,8)
Niter = 100
Nrep = 1000

params = {
    "p_bdy_qubit" : 0.1,
    "p_bulk_stab" : 0.01,
    "p_bdy_stab" : 0.1,
    "bandwidth" : 4,
    "code_dist" : 4,
    "dx" : 4,
    "repetitions": 4,
    "Niter" : Niter,
    "gen_coh_ratio" : gen_coh_ratio_list,
    "loss_prob" : 0.1
}


for i_d, d in enumerate(dist_list):
    dx = 2 * d
    dy = d
    repetitions = d
    code_details_2d = surface_code(dx, dy)
    code_details_3d = surface_code_3d(dx, dy,repetitions)
    code_details = code_details_2d, code_details_3d
    params["code_dist"] = d
    params["dx"] = dx
    params["repetitions"] = repetitions

    logical_err_rate = joint_meas_w_loss(code_details, params)
    print(logical_err_rate)
    # def runner(i_rep):
    #     tic = time.time()
    #     logical_err_rate = joint_meas_w_loss(code_details, params)
    #     out_dir = "results/joint_meas_loss/"
    #     fname = f"d_{d}_bw_{bandwidth}_r_{i_rep}.npz"
    #     toc = time.time()
    #     print(f"{fname}, elapsed time {toc-tic:.2f} sec")
    #     np.savez(out_dir + fname, gen_coh_ratio_list, logical_err_rate)
    # # runner(1)
    # results = Parallel(n_jobs=num_cores)(delayed(runner)(i_rep+364) for i_rep in range(Nrep))

print("Job finished!")
