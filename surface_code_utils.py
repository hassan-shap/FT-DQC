import numpy as np
import random
import time
from pymatching import Matching

def surface_code(dx, dy):
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

    return Smat, logicals, bdy, qubit_idx


def surface_code_3d(dx, dy,repetitions):

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

    qubit_bdy = [qubit_idx_spacetime[f"({dx/2+0.5}, {i_y}, {i_t})"] for i_y in range(dy) for i_t in range(repetitions)]
    S_bdy = [qubit_idx_spacetime[f"({dx//2}, {i_y}, {i_t+0.5})"] for i_y in range(dy) for i_t in range(repetitions-1)]
    bdy = qubit_bdy, S_bdy

    qubit_idx_mat = qubit_idx_spacetime, qubit_coords_spacetime, ancilla_idx

    # return Smat_spacetime, qubit_idx_mat
    return Smat_spacetime, bdy, qubit_idx_mat

harmonic = lambda k: (1/np.arange(1,k+1)).sum()

def perf_sim(d, bandwidth, code_details, gen_coh_ratio_list, p_err, Niter= 3000):
    
    p_bdy_qubit = p_err["bdy_qubit"]
    p_bulk_s = p_err["bulk_stab"]
    p_bdy_s = p_err["bdy_stab"]

    num_errors = np.zeros(len(gen_coh_ratio_list))
    Smat, logicals, bdy = code_details
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