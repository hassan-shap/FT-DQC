import numpy as np
import random
import time
from pymatching import Matching
import networkx as nx
import scipy as sp

def surface_code_star(dx, dy):
    # dx = 2 * d 
    # dy = d 

    qubit_coords = []
    for i_y in range(dy):
        for i_x in range(dx):
            if i_x ==0 or i_y == dy - 1:
                qubit_coords.append( ( i_x + 0.5, i_y ) )
            else:
                qubit_coords.append( ( i_x + 0.5, i_y ) )
                qubit_coords.append( ( i_x, i_y + 0.5 ) )

    
    qubit_idx = {}
    for idx, coords in enumerate(qubit_coords):
        qubit_idx[f"{coords}"] = idx

    # plotter(qubit_coords, color="C0")

    num_q = len(qubit_coords)
    num_s = (dx-1) * dy
    s_idx = 0
    Smat = np.zeros(( num_s, num_q))
    Smat_idx = {}
    for i_y in range(dy):
        for i_x in range(1,dx):
            Smat_idx[f"({i_x}, {i_y})"] = s_idx
            if i_y == 0 :
                qvec = [(i_x + 0.5, i_y), (i_x - 0.5, i_y), (i_x, i_y + 0.5)]
            elif i_y == dy - 1 :
                qvec = [(i_x + 0.5, i_y), (i_x - 0.5, i_y), (i_x, i_y - 0.5)]
            else:
                qvec = [(i_x + 0.5, i_y), (i_x - 0.5, i_y), (i_x, i_y - 0.5), (i_x, i_y + 0.5)]
            for q in qvec:
                Smat[s_idx, qubit_idx[f"{q}"]] = 1
            s_idx += 1

    qubit_bdy = [qubit_idx[f"({dx/2+0.5}, {i_y})"] for i_y in range(dy)]
    S_bdy = [Smat_idx[f"({dx//2}, {i_y})"] for i_y in range(dy)]
    bdy = qubit_bdy, S_bdy

    qubit_logical = [qubit_idx[f"({0.5}, {i_y})"] for i_y in range(dy)]
    logicals = np.zeros(num_q, dtype = np.uint8)
    logicals[qubit_logical] = 1

    # qubit_maps = qubit_idx, qubit_coords
    qubit_idx_mat = qubit_idx, Smat_idx
    return Smat, logicals, bdy, qubit_idx_mat

def surface_code_plaquette(dx, dy):

    qubit_coords = []
    for i_y in range(dy):
        for i_x in range(dx):
            if i_x ==0 or i_y == dy - 1:
                qubit_coords.append( ( i_x + 0.5, i_y ) )
            else:
                qubit_coords.append( ( i_x + 0.5, i_y ) )
                qubit_coords.append( ( i_x, i_y + 0.5 ) )
    
    qubit_idx = {}
    for idx, coords in enumerate(qubit_coords):
        qubit_idx[f"{coords}"] = idx

    num_q = len(qubit_coords)
    num_s = dx * (dy-1)
    s_idx = 0
    Smat = np.zeros(( num_s, num_q))
    Smat_idx = {}
    for i_y in range(dy-1):
        for i_x in range(1,1+dx):
            p_x = i_x-0.5
            p_y = i_y+0.5
            Smat_idx[f"({p_x}, {p_y})"] = s_idx
            qvec =[]
            if i_x == 1 :
                qvec = [(int(p_x + 0.5), p_y), (p_x, int(p_y -0.5)), (p_x, int(p_y + 0.5))]
            elif i_x == dx :
                qvec = [(int(p_x - 0.5), p_y), (p_x, int(p_y -0.5)), (p_x, int(p_y + 0.5))]
            else:
                qvec = [(int(p_x - 0.5), p_y), (int(p_x + 0.5), p_y), (p_x, int(p_y -0.5)), (p_x, int(p_y + 0.5))]
            for q in qvec:
                Smat[s_idx, qubit_idx[f"{q}"]] = 1
            s_idx += 1

    qubit_bdy = [qubit_idx[f"({dx//2}, {i_y+0.5})"] for i_y in range(dy-1)]
    S_bdy = [Smat_idx[f"({dx//2+0.5}, {i_y+0.5})"] for i_y in range(dy-1)]
    bdy = qubit_bdy, S_bdy

    qubit_logical = [qubit_idx[f"({i_x+0.5}, {0})"] for i_x in range(dx)]
    logicals = np.zeros(num_q)
    logicals[qubit_logical] = 1

    return Smat, logicals, bdy, qubit_idx


def surface_code_3d_star(dx, dy,repetitions):

    # dx = 2*d
    # dy = d 
    # repetitions = d

    qubit_coords_spacetime = []
    for i_t in range(repetitions):
        for i_y in range(dy):
            for i_x in range(dx):
                if i_t < repetitions - 1:
                    if i_x ==0:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                    elif i_y == dy - 1:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                        qubit_coords_spacetime.append( ( i_x, i_y, i_t + 0.5 ) )
                    else:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                        qubit_coords_spacetime.append( ( i_x, i_y + 0.5, i_t ) )
                        qubit_coords_spacetime.append( ( i_x, i_y, i_t + 0.5 ) )
                else:
                    if i_x ==0:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                    elif i_y == dy - 1:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                    else:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                        qubit_coords_spacetime.append( ( i_x, i_y + 0.5, i_t ) )

    
    qubit_idx_spacetime = {}
    # qubit_idx_to_coords_spacetime = {}
    for idx, coords in enumerate(qubit_coords_spacetime):
        qubit_idx_spacetime[f"{coords}"] = idx
        # qubit_idx_to_coords_spacetime[idx] = coords

    num_q_spacetime = len(qubit_coords_spacetime)
    num_s_spacetime = (dx-1) * dy * repetitions
    s_idx = 0
    Smat_spacetime = np.zeros(( num_s_spacetime, num_q_spacetime))
    Smat_idx_spacetime = {}
    for i_t in range(repetitions):
        for i_y in range(dy):
            for i_x in range(1,dx):
                Smat_idx_spacetime[f"({i_x}, {i_y}, {i_t})"] = s_idx
                if i_t == 0:
                    if i_y == 0 :
                        qvec = [(i_x + 0.5, i_y, i_t), (i_x - 0.5, i_y, i_t), (i_x, i_y + 0.5, i_t),
                                (i_x, i_y, i_t + 0.5)
                                ]
                    elif i_y == dy - 1 :
                        qvec = [(i_x + 0.5, i_y, i_t), (i_x - 0.5, i_y, i_t), (i_x, i_y - 0.5, i_t),
                                (i_x, i_y, i_t + 0.5)
                                ]
                    else:
                        qvec = [(i_x + 0.5, i_y, i_t), (i_x - 0.5, i_y, i_t), (i_x, i_y + 0.5, i_t),
                                (i_x, i_y - 0.5, i_t), (i_x, i_y, i_t + 0.5)
                                ]
                elif i_t == repetitions - 1:
                    if i_y == 0 :
                        qvec = [(i_x + 0.5, i_y, i_t), (i_x - 0.5, i_y, i_t), (i_x, i_y + 0.5, i_t),
                                (i_x, i_y, i_t - 0.5)
                                ]
                    elif i_y == dy - 1 :
                        qvec = [(i_x + 0.5, i_y, i_t), (i_x - 0.5, i_y, i_t), (i_x, i_y - 0.5, i_t),
                                (i_x, i_y, i_t - 0.5)
                                ]
                    else:
                        qvec = [(i_x + 0.5, i_y, i_t), (i_x - 0.5, i_y, i_t), (i_x, i_y + 0.5, i_t),
                                (i_x, i_y - 0.5, i_t), (i_x, i_y, i_t - 0.5)
                                ]
                else:
                    if i_y == 0 :
                        qvec = [(i_x + 0.5, i_y, i_t), (i_x - 0.5, i_y, i_t), (i_x, i_y + 0.5, i_t),
                                (i_x, i_y, i_t + 0.5), (i_x, i_y, i_t - 0.5)
                                ]
                    elif i_y == dy - 1 :
                        qvec = [(i_x + 0.5, i_y, i_t), (i_x - 0.5, i_y, i_t), (i_x, i_y - 0.5, i_t),
                                (i_x, i_y, i_t + 0.5), (i_x, i_y, i_t - 0.5)
                                ]
                    else:
                        qvec = [(i_x + 0.5, i_y, i_t), (i_x - 0.5, i_y, i_t), (i_x, i_y + 0.5, i_t),
                                (i_x, i_y - 0.5, i_t), (i_x, i_y, i_t + 0.5), (i_x, i_y, i_t - 0.5)
                                ]
                for q in qvec:
                    Smat_spacetime[s_idx, qubit_idx_spacetime[f"{q}"]] = 1
                s_idx += 1

    ancilla_idx = [qubit_idx_spacetime[f"({i_x}, {i_y}, {i_t+0.5})"] for i_x in range(1,dx) for i_y in range(dy) for i_t in range(repetitions-1)]

    qubit_bdy = sorted([qubit_idx_spacetime[f"({dx/2+0.5}, {i_y}, {i_t})"] for i_y in range(dy) for i_t in range(repetitions)])
    S_bdy = sorted([qubit_idx_spacetime[f"({dx//2}, {i_y}, {i_t+0.5})"] for i_y in range(dy) for i_t in range(repetitions-1)])
    bdy = qubit_bdy, S_bdy

    qubit_idx_mat = qubit_idx_spacetime, qubit_coords_spacetime, ancilla_idx

    # return Smat_spacetime, qubit_idx_mat
    return Smat_spacetime, bdy, qubit_idx_mat

def surface_code_3d_plaquette(dx, dy,repetitions):

    qubit_coords_spacetime = []
    for i_t in range(repetitions):
        if i_t < repetitions - 1:
            for i_y in range(dy):
                for i_x in range(dx):
                    if i_x ==0:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                    elif i_y == dy - 1:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                    else:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                        qubit_coords_spacetime.append( ( i_x, i_y + 0.5, i_t ) )

            for i_y in range(dy-1):
                for i_x in range(1,1+dx):
                    p_x = i_x-0.5
                    p_y = i_y+0.5
                    qubit_coords_spacetime.append( ( p_x, p_y, i_t + 0.5) )
        else:
            for i_y in range(dy):
                for i_x in range(dx):
                    if i_x ==0:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                    elif i_y == dy - 1:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                    else:
                        qubit_coords_spacetime.append( ( i_x + 0.5, i_y, i_t ) )
                        qubit_coords_spacetime.append( ( i_x, i_y + 0.5, i_t ) )
    
    qubit_idx_spacetime = {}
    ancilla_idx = []
    for idx, coords in enumerate(qubit_coords_spacetime):
        qubit_idx_spacetime[f"{coords}"] = idx
        if coords[2] % 1 > 0:
            ancilla_idx.append(idx)

    num_q_spacetime = len(qubit_coords_spacetime)
    num_s_spacetime = dx * (dy-1) * repetitions
    s_idx = 0
    Smat_spacetime = np.zeros(( num_s_spacetime, num_q_spacetime))
    Smat_idx_spacetime = {}
    for i_t in range(repetitions):
        for i_y in range(dy-1):
            for i_x in range(1,1+dx):
                p_x = i_x-0.5
                p_y = i_y+0.5
                Smat_idx_spacetime[f"({p_x}, {p_y}, {i_t})"] = s_idx
                qvec =[]
                if i_t == 0:
                    if i_x == 1 :
                        qvec = [(int(p_x + 0.5), p_y, i_t), (p_x, int(p_y -0.5), i_t), (p_x, int(p_y + 0.5), i_t),
                                (p_x, p_y, i_t + 0.5)
                                ]
                    elif i_x == dx :
                        qvec = [(int(p_x - 0.5), p_y, i_t), (p_x, int(p_y -0.5), i_t), (p_x, int(p_y + 0.5), i_t),
                                (p_x, p_y, i_t + 0.5)
                                ]
                    else:
                        qvec = [(int(p_x - 0.5), p_y, i_t), (int(p_x + 0.5), p_y, i_t), (p_x, int(p_y -0.5), i_t),
                                (p_x, int(p_y + 0.5), i_t), (p_x, p_y, i_t + 0.5)
                                ]
                elif i_t == repetitions - 1:
                    if i_x == 1 :
                        qvec = [(int(p_x + 0.5), p_y, i_t), (p_x, int(p_y -0.5), i_t), (p_x, int(p_y + 0.5), i_t),
                                (p_x, p_y, i_t - 0.5)
                                ]
                    elif i_x == dx :
                        qvec = [(int(p_x - 0.5), p_y, i_t), (p_x, int(p_y -0.5), i_t), (p_x, int(p_y + 0.5), i_t),
                                (p_x, p_y, i_t - 0.5)
                                ]
                    else:
                        qvec = [(int(p_x - 0.5), p_y, i_t), (int(p_x + 0.5), p_y, i_t), (p_x, int(p_y -0.5), i_t),
                                (p_x, int(p_y + 0.5), i_t), (p_x, p_y, i_t - 0.5)
                                ]
                else:
                    if i_x == 1 :
                        qvec = [(int(p_x + 0.5), p_y, i_t), (p_x, int(p_y -0.5), i_t), (p_x, int(p_y + 0.5), i_t),
                                (p_x, p_y, i_t + 0.5), (p_x, p_y, i_t - 0.5)
                                ]
                    elif i_x == dx :
                        qvec = [(int(p_x - 0.5), p_y, i_t), (p_x, int(p_y -0.5), i_t), (p_x, int(p_y + 0.5), i_t),
                                (p_x, p_y, i_t + 0.5), (p_x, p_y, i_t - 0.5)
                                ]
                    else:
                        qvec = [(int(p_x - 0.5), p_y, i_t), (int(p_x + 0.5), p_y, i_t), (p_x, int(p_y -0.5), i_t),
                                (p_x, int(p_y + 0.5), i_t), (p_x, p_y, i_t + 0.5), (p_x, p_y, i_t - 0.5)
                                ]
                for q in qvec:
                    Smat_spacetime[s_idx, qubit_idx_spacetime[f"{q}"]] = 1
                s_idx += 1

    # ancilla_idx = [qubit_idx_spacetime[f"({i_x}, {i_y}, {i_t+0.5})"] for i_x in range(1,dx) for i_y in range(dy) for i_t in range(repetitions-1)]

    qubit_bdy = sorted([qubit_idx_spacetime[f"({dx//2}, {i_y+0.5}, {i_t})"] for i_y in range(dy-1) for i_t in range(repetitions)])
    S_bdy = sorted([qubit_idx_spacetime[f"({dx/2+0.5}, {i_y+0.5}, {i_t+0.5})"] for i_y in range(dy-1) for i_t in range(repetitions-1)])
    bdy = qubit_bdy, S_bdy

    qubit_idx_mat = qubit_idx_spacetime, qubit_coords_spacetime, ancilla_idx
    # logicals =  logical_vertical, logical_horizontal

    return Smat_spacetime, bdy, qubit_idx_mat#, logicals

harmonic = lambda k: (1/np.arange(1,k+1)).sum()

def perf_sim(d, bandwidth, code_details, gen_coh_ratio_list, p_err, Niter= 3000):
    
    p_bdy_qubit = p_err["bdy_qubit"]
    p_bulk_s = p_err["bulk_stab"]
    p_bdy_s = p_err["bdy_stab"]

    num_errors = np.zeros(len(gen_coh_ratio_list))
    Smat, logicals, bdy, _ = code_details
    num_s, num_q = Smat.shape
    qubit_bdy, S_bdy = bdy
    p_s_list = np.ones(num_s) * p_bulk_s
    p_s_list[S_bdy] = p_bdy_s

    repetitions = d

    for i_t, ratio in enumerate(gen_coh_ratio_list):
        # print(i_t)
        
        if d > bandwidth:
            idle_time = ((d // bandwidth) * harmonic(bandwidth)  + harmonic (d % bandwidth)  )* ratio
        else:
            idle_time = harmonic(d)* ratio
        p_bulk_qubit = 1 - np.exp(- idle_time )
        p_qubit_list = np.ones(num_q) * p_bulk_qubit
        p_qubit_list[qubit_bdy] = (p_bulk_qubit+ p_bdy_qubit)

        matching = Matching(Smat, spacelike_weights=np.log((1-p_qubit_list)/p_qubit_list),
                        repetitions=repetitions, timelike_weights=np.log((1-p_s_list)/p_s_list), faults_matrix=logicals)
        for _ in range(Niter):
            noise_new = (np.random.rand(num_q, repetitions) < np.reshape(p_qubit_list,(p_qubit_list.shape[0],1))).astype(np.uint8)
            noise_cumulative = (np.cumsum(noise_new, 1) % 2).astype(np.uint8)
            noise_total = noise_cumulative[:,-1]
            syndrome = Smat@noise_cumulative % 2
            syndrome_error = (np.random.rand(num_s, repetitions) < np.reshape(p_s_list,(p_s_list.shape[0],1)) ).astype(np.uint8)
            syndrome_error[:,-1] = 0 # Perfect measurements in last round to ensure even parity
            noisy_syndrome = (syndrome + syndrome_error) % 2
            # Convert to difference syndrome
            noisy_syndrome[:,1:] = (noisy_syndrome[:,1:] - noisy_syndrome[:,0:-1]) % 2
            # recovery_chain = matching.decode(noisy_syndrome)
            # print(num_q, num_s, recovery_chain.shape)
            predicted_logicals_flipped = matching.decode(noisy_syndrome)[0]
            actual_logicals_flipped = noise_total@logicals.T % 2
            num_errors[i_t] += not (predicted_logicals_flipped == actual_logicals_flipped)
    num_errors /= Niter
    # print(num_errors)
    # print(predicted_logicals_flipped, actual_logicals_flipped)
    return num_errors

def joint_meas_w_loss(code_details, params):

    p_bdy_qubit = params["p_bdy_qubit"]
    p_bulk_s = params["p_bulk_stab"]
    p_bdy_s = params["p_bdy_stab"]
    bandwidth = params["bandwidth"]
    repetitions = params["repetitions"]
    Niter = params["Niter"]
    d = params["code_dist"]
    dx = params["dx"]
    gen_coh_ratio_list = params["gen_coh_ratio"]
    num_errors = np.zeros(len(gen_coh_ratio_list))

    code_details_2d, code_details_3d = code_details
    Smat, logicals, bdy, qubit_idx = code_details_2d
    num_s, num_q = Smat.shape
    data_bdy, ancilla_bdy = bdy

    Smat_spacetime, bdy_spacetime, qubit_idx_mat = code_details_3d
    data_bdy_st, ancilla_bdy_st = bdy_spacetime
    qubit_idx_spacetime, qubit_coords_spacetime, ancilla_idx_st = qubit_idx_mat
    num_q_spacetime = Smat_spacetime.shape[1]
    data_q_idx = list(set(range(num_q_spacetime))-set(ancilla_idx_st))

    q_idx_mat = np.zeros( (repetitions, num_q), dtype= np.uint32)
    for idx in range(num_q_spacetime):
        x, y, t = qubit_coords_spacetime[idx]
        if t % 1 == 0:
            q_idx_mat[t, qubit_idx[f"({x}, {y})"]] = idx

    # ############## error syndrome calculations ######

    p_s_list = np.ones(num_s) * p_bulk_s
    p_s_list[ancilla_bdy] = p_bdy_s

    for i_t, ratio in enumerate(gen_coh_ratio_list):
        # print(i_t)
        
        if d > bandwidth:
            idle_time = ((d // bandwidth) * harmonic(bandwidth)  + harmonic (d % bandwidth)  )* ratio
        else:
            idle_time = harmonic(d)* ratio
        p_bulk_qubit = 1 - np.exp(- idle_time )
        p_qubit_list = np.ones(num_q) * p_bulk_qubit
        p_qubit_list[data_bdy] = (p_bulk_qubit+ p_bdy_qubit)

        p_qubit_st = np.zeros(num_q_spacetime) 
        p_qubit_st[data_q_idx] = p_bulk_qubit 
        p_qubit_st[data_bdy_st] = (p_bulk_qubit+ p_bdy_qubit) 
        p_qubit_st[ancilla_idx_st] = p_bulk_s
        p_qubit_st[ancilla_bdy_st] = p_bdy_s

        for _ in range(Niter):

            prob_l = params["loss_prob"]
            error_loss = np.random.rand(num_s , repetitions)
            error_loss[:,-1] = 1
            loss_inds_2d = np.argwhere(error_loss < prob_l)
            loss_inds = np.array([qubit_idx_spacetime[f"({1+ (coords[0]% (dx-1))}, {coords[0]//(dx-1)}, {coords[1]+0.5})"] for coords in loss_inds_2d])
            remain_inds = np.array(list( set(range(num_q_spacetime)) - set(loss_inds)))

            Sx_red, qubits_to_plot = compute_eff_Sx(Smat_spacetime,loss_inds,remain_inds)
            overlap = Sx_red.T@Sx_red
            ql_rep = find_repeated_idx(overlap,qubits_to_plot)

            remain_ancilla = np.intersect1d(remain_inds, ancilla_idx_st)
            remain_ancilla_bdy = np.intersect1d(remain_inds, ancilla_bdy_st)

            weights = np.zeros(num_q_spacetime)
            weights[remain_ancilla] = np.log((1-p_bulk_s)/p_bulk_s) 
            weights[remain_ancilla_bdy] = np.log((1-p_bdy_s)/p_bdy_s) 
            weights[data_q_idx] = np.log((1-p_bulk_qubit)/p_bulk_qubit) 
            weights[data_bdy_st] = np.log((1-p_bdy_qubit)/p_bdy_qubit)

            for links in ql_rep:
                p_eff = (1-np.prod(1-2*p_qubit_st[links]))/2
                weights[links] = np.log((1-p_eff)/p_eff)

            assert len(np.argwhere(weights>0))== len(qubits_to_plot)

            matching_st = Matching(Smat_spacetime,spacelike_weights=weights)

            noise_new = (np.random.rand(num_q, repetitions) < np.reshape(p_qubit_list,(p_qubit_list.shape[0],1))).astype(np.uint8)
            noise_cumulative = (np.cumsum(noise_new, 1) % 2).astype(np.uint8)
            noise_total = noise_cumulative[:,-1]
            syndrome = Smat@noise_cumulative % 2

            syndrome_error = (np.random.rand(num_s, repetitions) < np.reshape(p_s_list,(p_s_list.shape[0],1)) ).astype(np.uint8)
            syndrome_error[:,-1] = 0 # Perfect measurements in last round to ensure even parity
            noisy_syndrome = (syndrome + syndrome_error) % 2
            # # Convert to difference syndrome
            noisy_syndrome[:,1:] = (noisy_syndrome[:,1:] - noisy_syndrome[:,0:-1]) % 2

            recovery_st = matching_st.decode(noisy_syndrome)
            # recovery_tot = np.zeros(num_q)
            # for idx, q in enumerate(recovery_st):
            #     x, y, t = qubit_coords_spacetime[idx]
            #     if t % 1 == 0:
            #         recovery_tot[qubit_idx[f"({x}, {y})"]] += recovery_st[idx]
            # recovery_tot %= 2
            recovery_tot = np.zeros(num_q)
            for t in range(repetitions):
                recovery_tot += recovery_st[q_idx_mat[t,:]]
            recovery_tot %= 2

            predicted_logicals_flipped = recovery_tot@logicals.T % 2
            actual_logicals_flipped = noise_total@logicals.T % 2
            num_errors[i_t] += not (predicted_logicals_flipped == actual_logicals_flipped)

    num_errors /= Niter

    return num_errors

def perf_sim_no_loss(code_details, params):

    p_bdy_qubit = params["p_bdy_qubit"]
    p_bulk_s = params["p_bulk_stab"]
    p_bdy_s = params["p_bdy_stab"]
    repetitions = params["repetitions"]
    Niter = params["Niter"]
    d = params["code_dist"]
    dx = params["dx"]
    p_bulk_qubit_list = params["p_bulk_qubit"]
    num_errors = np.zeros(len(p_bulk_qubit_list))
    
    Smat, logicals, bdy, _ = code_details
    num_s, num_q = Smat.shape
    qubit_bdy, S_bdy = bdy
    p_s_list = np.ones(num_s) * p_bulk_s
    p_s_list[S_bdy] = p_bdy_s

    for i_p, p_bulk_qubit in enumerate(p_bulk_qubit_list):
        p_qubit_list = np.ones(num_q) * p_bulk_qubit
        p_qubit_list[qubit_bdy] = (p_bulk_qubit+ p_bdy_qubit)

        matching = Matching(Smat, spacelike_weights=np.log((1-p_qubit_list)/p_qubit_list),
                        repetitions=repetitions, timelike_weights=np.log((1-p_s_list)/p_s_list), faults_matrix=logicals)
        for _ in range(Niter):
            noise_new = (np.random.rand(num_q, repetitions) < np.reshape(p_qubit_list,(p_qubit_list.shape[0],1))).astype(np.uint8)
            noise_cumulative = (np.cumsum(noise_new, 1) % 2).astype(np.uint8)
            noise_total = noise_cumulative[:,-1]
            syndrome = Smat@noise_cumulative % 2
            syndrome_error = (np.random.rand(num_s, repetitions) < np.reshape(p_s_list,(p_s_list.shape[0],1)) ).astype(np.uint8)
            syndrome_error[:,-1] = 0 # Perfect measurements in last round to ensure even parity
            noisy_syndrome = (syndrome + syndrome_error) % 2
            # Convert to difference syndrome
            noisy_syndrome[:,1:] = (noisy_syndrome[:,1:] - noisy_syndrome[:,0:-1]) % 2
            # recovery_chain = matching.decode(noisy_syndrome)
            # print(num_q, num_s, recovery_chain.shape)
            predicted_logicals_flipped = matching.decode(noisy_syndrome)[0]
            actual_logicals_flipped = noise_total@logicals.T % 2
            num_errors[i_p] += not (predicted_logicals_flipped == actual_logicals_flipped)
    num_errors /= Niter
    # print(num_errors)
    # print(predicted_logicals_flipped, actual_logicals_flipped)
    return num_errors

def joint_meas_2d_phase_diagram(code_details, params):

    p_bulk_qubit_list = params["p_bulk_qubit"]
    num_errors = np.zeros(len(p_bulk_qubit_list))
    p_bdy_qubit = params["p_bdy_qubit"]
    p_bulk_s = params["p_bulk_stab"]
    p_bdy_s = params["p_bdy_stab"]
    prob_l = params["loss_prob"]
    repetitions = params["repetitions"]
    Niter = params["Niter"]
    dx = params["dx"]

    code_details_2d, code_details_3d = code_details
    Smat, logicals, bdy, qubit_idx = code_details_2d
    num_s, num_q = Smat.shape
    data_bdy, ancilla_bdy = bdy

    Smat_spacetime, bdy_spacetime, qubit_idx_mat = code_details_3d
    data_bdy_st, ancilla_bdy_st = bdy_spacetime
    qubit_idx_spacetime, qubit_coords_spacetime, ancilla_idx_st = qubit_idx_mat
    num_q_spacetime = Smat_spacetime.shape[1]
    data_q_idx = list(set(range(num_q_spacetime))-set(ancilla_idx_st))

    q_idx_mat = np.zeros( (repetitions, num_q), dtype= np.uint32)
    for idx in range(num_q_spacetime):
        x, y, t = qubit_coords_spacetime[idx]
        if t % 1 == 0:
            q_idx_mat[t, qubit_idx[f"({x}, {y})"]] = idx

    # ############## error syndrome calculations ######

    p_s_list = np.ones(num_s) * p_bulk_s
    p_s_list[ancilla_bdy] = p_bdy_s

    for i_p, p_bulk_qubit in enumerate(p_bulk_qubit_list):
        # print(i_t)
        
        p_qubit_list = np.ones(num_q) * p_bulk_qubit
        p_qubit_list[data_bdy] = (p_bulk_qubit+ p_bdy_qubit)

        p_qubit_st = np.zeros(num_q_spacetime) 
        p_qubit_st[data_q_idx] = p_bulk_qubit 
        p_qubit_st[data_bdy_st] = (p_bulk_qubit+ p_bdy_qubit) 
        p_qubit_st[ancilla_idx_st] = p_bulk_s
        p_qubit_st[ancilla_bdy_st] = p_bdy_s

        for _ in range(Niter):

            # error_loss = np.random.rand(len(data_bdy_st) , repetitions)
            # error_loss[:,-1] = 1
            # loss_inds_2d = data_bdy_st[np.argwhere(error_loss < prob_l)]
            # loss_inds = np.array([qubit_idx_spacetime[f"({1+ (coords[0]% (dx-1))}, {coords[0]//(dx-1)}, {coords[1]+0.5})"] for coords in loss_inds_2d])
            err_loss = np.argwhere(np.random.rand(len(ancilla_bdy_st))<prob_l)[:,0]
            loss_inds = [ancilla_bdy_st[idx] for idx in err_loss]
            remain_inds = np.array(list( set(range(num_q_spacetime)) - set(loss_inds)))

            # tic = time.time()

            Sx_red, qubits_to_plot = compute_eff_Sx(Smat_spacetime,loss_inds,remain_inds)
            # # toc = time.time()
            # print(f"s_eff takes {toc-tic:.2f} sec.")
            # # Sx_red = sp.sparse.coo_array(Sx_red)
            # tic = time.time()
            overlap = Sx_red.T@Sx_red
            # # toc = time.time()
            # # print(f"rep idx takes {toc-tic:.2f} sec.")
            # # tic = time.time()
            # G = bipartite.from_biadjacency_matrix(sp.sparse.csr_array(Sx_red))
            # X, Y = bipartite.sets(G)
            # nodes = [ n for n in X if G.degree[n]>6]
            # overlap = Sx_red[nodes,:].T@Sx_red[nodes,:]
            # # toc = time.time()
            # # print(f"rep idx eff takes {toc-tic:.2f} sec.")
            ql_rep = find_repeated_idx(overlap,qubits_to_plot)
            # tic = time.time()
            # ql_rep = find_repeated_idx_eff(Sx_red,qubits_to_plot)
            # toc = time.time()
            # print(f"rep idx eff takes {toc-tic:.2f} sec.")
            
            remain_ancilla = np.intersect1d(remain_inds, ancilla_idx_st)
            remain_ancilla_bdy = np.intersect1d(remain_inds, ancilla_bdy_st)

            weights = np.zeros(num_q_spacetime)
            weights[remain_ancilla] = np.log((1-p_bulk_s)/p_bulk_s) 
            weights[remain_ancilla_bdy] = np.log((1-p_bdy_s)/p_bdy_s) 
            weights[data_q_idx] = np.log((1-p_bulk_qubit)/p_bulk_qubit) 
            weights[data_bdy_st] = np.log((1-p_bdy_qubit)/p_bdy_qubit)

            for links in ql_rep:
                p_eff = (1-np.prod(1-2*p_qubit_st[links]))/2
                weights[links] = np.log((1-p_eff)/p_eff)

            assert len(np.argwhere(weights>0))== len(qubits_to_plot)

            matching_st = Matching(Smat_spacetime,spacelike_weights=weights)

            noise_new = (np.random.rand(num_q, repetitions) < np.reshape(p_qubit_list,(p_qubit_list.shape[0],1))).astype(np.uint8)
            noise_cumulative = (np.cumsum(noise_new, 1) % 2).astype(np.uint8)
            noise_total = noise_cumulative[:,-1]
            syndrome = Smat@noise_cumulative % 2

            syndrome_error = (np.random.rand(num_s, repetitions) < np.reshape(p_s_list,(p_s_list.shape[0],1)) ).astype(np.uint8)
            syndrome_error[:,-1] = 0 # Perfect measurements in last round to ensure even parity
            noisy_syndrome = (syndrome + syndrome_error) % 2
            ## Convert to difference syndrome
            noisy_syndrome[:,1:] = (noisy_syndrome[:,1:] - noisy_syndrome[:,0:-1]) % 2

            recovery_st = matching_st.decode(noisy_syndrome)
            # recovery_tot = np.zeros(num_q)
            # for idx, q in enumerate(recovery_st):
            #     x, y, t = qubit_coords_spacetime[idx]
            #     if t % 1 == 0:
            #         recovery_tot[qubit_idx[f"({x}, {y})"]] += recovery_st[idx]
            # recovery_tot %= 2

            recovery_tot = np.zeros(num_q)
            for t in range(repetitions):
                recovery_tot += recovery_st[q_idx_mat[t,:]]

            predicted_logicals_flipped = recovery_tot@logicals.T % 2
            actual_logicals_flipped = noise_total@logicals.T % 2
            num_errors[i_p] += not (predicted_logicals_flipped == actual_logicals_flipped)

    num_errors /= Niter

    return num_errors


##################
def compute_eff_Sx(Sx,loss_inds,remain_inds):
    G_loss = nx.Graph()
    for loss_index in loss_inds:
        s_inds = np.argwhere(Sx[:,loss_index]>0)[:,0]
        # if len(s_inds)>1:
        G_loss.add_edge(s_inds[0],s_inds[1])
        # else:
        #     G_loss.add_node(s_inds[0])

    components = [G_loss.subgraph(c).copy() for c in nx.connected_components(G_loss)]
    lost_vs = []
    for i_c, c in enumerate(components):
        lost_vs += c.nodes()

    num_s = Sx.shape[0]
    remain_vs = list(set(range(num_s)) - set(lost_vs))
    num_stab = len(components)+len(remain_vs)
    Sx_red2 = np.zeros((num_stab,len(remain_inds)))
    Sx_red2[len(components):,:] = Sx[np.ix_(remain_vs,remain_inds)]
    for i_c, c in enumerate(components):
        Sx_red2[i_c,:] = np.sum(Sx[np.ix_(c.nodes(),remain_inds)],axis = 0)%2

    keep_cols = np.argwhere(np.sum(Sx_red2,axis=0)>0)[:,0]
    Sx_red2 = Sx_red2[:,keep_cols]
    qubits_to_plot = remain_inds[keep_cols]

    return Sx_red2, qubits_to_plot

def find_repeated_idx(overlap,qubits_to_plot):
    inds = np.argwhere(overlap>1)
    if len(inds)>0:
        rep_edges = []
        for i_v in inds:
            if i_v[1]>i_v[0]:
                if not (i_v[0] in rep_edges):
                    rep_edges.append(i_v[0])
                if not (i_v[1] in rep_edges):
                    rep_edges.append(i_v[1])

        rep_edges = np.sort(rep_edges).astype(int)

        repeated_inds = []
        counter = 0
        i = 0 
        overlap2 = overlap[np.ix_(rep_edges,rep_edges)]
        inds_to_keep2 = list(range(len(rep_edges)))
        while counter < len(rep_edges):
            edge = inds_to_keep2[i]
            ovlp_inds = np.argwhere(overlap2[edge,inds_to_keep2[i+1:]]==2)
            nl_i = len(ovlp_inds)+1

            qlist = qubits_to_plot[rep_edges[np.ix_([ inds_to_keep2[k] for k in i+1+ovlp_inds[:,0]])]]
            for j in ovlp_inds[::-1,0]:
                inds_to_keep2.remove(inds_to_keep2[i+1+j])
            repeated_inds.append(list(np.concatenate(([qubits_to_plot[rep_edges[edge]]],qlist))))
            counter += nl_i
            i += 1

    else:
        repeated_inds = []

    return repeated_inds
# ##################

def find_repeated_idx_eff(overlap,qubits_to_plot):
    num_qubits = len(qubits_to_plot)
    overlap[range(num_qubits),range(num_qubits)] = 0
    inds = np.argwhere(overlap>1)
    if len(inds)>1:
        repeated_inds = []
        for i_v in inds:
            idx2 = list(np.argwhere(overlap[i_v[0],:]>1)[:,0])
            repeated_inds.append(list(qubits_to_plot[ [i_v[0]] + idx2]))
            
    else:
        repeated_inds = []
    
    return repeated_inds

from networkx.algorithms import bipartite

def find_repeated_idx_eff_2(Sx_red,qubits_to_plot):

    G = bipartite.from_biadjacency_matrix(sp.sparse.csr_array(Sx_red))
    X, Y = bipartite.sets(G)
    nodes = [ n for n in X if G.degree[n]>6]
    overlap = Sx_red[nodes,:].T@Sx_red[nodes,:]
    num_qubits = len(qubits_to_plot)
    overlap[range(num_qubits),range(num_qubits)] = 0
    inds = np.argwhere(overlap>1)
    if len(inds)>1:
        repeated_inds = []
        for i_v in inds:
            idx2 = list(np.argwhere(overlap[i_v[0],:]>1)[:,0])
            repeated_inds.append(list(qubits_to_plot[ [i_v[0]] + idx2]))
            
    else:
        repeated_inds = []
    
    return repeated_inds

