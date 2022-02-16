import numpy as np
import matplotlib.pyplot as plt
from math import pi
import sys
import os.path

DATA_DIR = os.path.dirname(__file__)+'/../Data/'
OUT_DIR = './figs/'
theta_pi_frac = 1.0
theta = theta_pi_frac*pi
nsets = 100
nstart = 10
nproc = 2
injmod = 1
num_steps_log = 200
E_inj_exp = 10
#stepexp = 2.1
#stepsize = 0

flags = ['--theta-pi-fr', '--tmax', '--stepexp']

if __name__ == "__main__":
    if (len(sys.argv)-1)%2 == 0:
        args = sys.argv[1:]
        for i in range(0, len(args), 2):
            flag = args[i]
            val = args[i+1]

            if flag in flags:
                for j, f in enumerate(flags):
                    if flag == f:
                        if j == 0:
                            try:
                                theta_pi_frac = float(val)
                                theta = theta_pi_frac * pi
                                print(f"Using provided theta_pi_frac{theta_pi_frac}")
                            except:
                                print(f"Invalid argument, using default theta_pi_frac{theta_pi_frac}    ")
                        elif j == 1:
                            try:
                                t_max = int(val)
                                print(f"Using provided t_max{t_max}")
                            except:
                                print(f"Invalid argument, using default t_max{t_max}")
                        elif j == 2:
                            try:
                                stepexp = float(val)
                                print(f"Using provided stepexp{stepexp}")
                            except:
                                print(f"Invalid argument, using default stepexp{stepexp}")
            else:
                print(f"Bad flag {flag} provided, no new parameter value set")
    else:
        print("Odd number of flags+parameters, using default values")

    print("########################")
    print("theta_pi_frac and theta")
    print(theta_pi_frac)
    print(theta)
    print("########################")
    print()

    phase_space_fname = DATA_DIR
    phase_space_fname += (
            f"phase_space_acc_injmod{injmod}_"
            f"theta-max{theta:.3f}_"
            f"nsets{nsets}_"
            f"nstart{nstart}_"
            f"E_inj_exp{E_inj_exp:.3f}_"
            f"nproc{nproc}"
        )
    phase_space_file = open(phase_space_fname, 'rb')
    phase_space_data = np.fromfile(phase_space_file)

    n_log = 200
    data_array_size = nstart*2*n_log
    set_size = data_array_size * nsets
    total_data_size = set_size * nproc
    print('##################')
    print("ERROR: ")
    print("data array size: ", data_array_size) # Correct
    print("set size: ", set_size) # Should be correct
    print("Expected total size: ", total_data_size)
    print("Actual size: ", phase_space_data.size)
    print('##################')
    raise Exception



    d_sampled = np.zeros(n_samples)
    p_sampled = np.zeros(n_samples)

    for sample in range(n_samples):
        d_samples = np.zeros(nsets*nstart*nproc)
        p_samples = np.zeros(nsets*nstart*nproc)
        for seti in range(nsets):
            d_samples_seti = np.zeros(nstart*nproc)
            p_samples_seti = np.zeros(nstart*nproc)
            for proc in range(nproc):
                data = phase_space_data[proc, seti, sample]
                d, p = data.T
                d_samples_seti[proc*nstart:(proc+1)*nstart] = d
                p_samples_seti[proc*nstart:(proc+1)*nstart] = p
                if proc == 0 and seti == 0:
                    t_sampled[sample] = t[0]

    phase_space_fname_iso = DATA_DIR
    phase_space_fname_iso += (
            f"phase_space_rwk_tmax{t_max:.3f}_"
            f"theta-max{np.pi:.3f}_"
            f"nsets{nsets}_"
            f"nstart{nstart}_"
            f"E_inj_exp{E_inj_exp:.3f}_"
            f"nproc{nproc}"
        )

    avg_drifts_sampled = np.zeros(n_samples)
    t_sampled = np.zeros(n_samples)

    for sample in range(n_samples):
        x_samples = np.zeros(nsets*nstart*nproc)
        y_samples = np.zeros(nsets*nstart*nproc)
        z_samples = np.zeros(nsets*nstart*nproc)
        for seti in range(nsets):
            x_samples_seti = np.zeros(nstart*nproc)
            y_samples_seti = np.zeros(nstart*nproc)
            z_samples_seti = np.zeros(nstart*nproc)
            for proc in range(nproc):
                data = phase_space_data[proc, seti, sample]
                x, y, z, t = data.T
                x_samples_seti[proc*nstart:(proc+1)*nstart] = x
                y_samples_seti[proc*nstart:(proc+1)*nstart] = y
                z_samples_seti[proc*nstart:(proc+1)*nstart] = z
                if proc == 0 and seti == 0:
                    t_sampled[sample] = t[0]
            x_samples[seti*len(x_samples_seti):(seti+1)*len(x_samples_seti)] = x_samples_seti
            y_samples[seti*len(y_samples_seti):(seti+1)*len(y_samples_seti)] = y_samples_seti
            z_samples[seti*len(z_samples_seti):(seti+1)*len(z_samples_seti)] = z_samples_seti
        avg_drifts_sampled[sample] = np.average(x_samples**2 + y_samples**2 + z_samples**2)

    # isotropic data
    phase_space_fname_iso = DATA_DIR
    phase_space_fname_iso += (
            f"phase_space_rwk_tmax{t_max:.3f}_"
            f"theta-max{np.pi:.3f}_"
            f"nsets{nsets}_"
            f"nstart{nstart}_"
            f"E_inj_exp{E_inj_exp:.3f}_"
            f"nproc{nproc}"
        )
    phase_space_file_iso = open(phase_space_fname_iso, 'rb')
    phase_space_data_iso = np.fromfile(phase_space_file_iso, dtype='float64')
    n_samples_iso = int(len(phase_space_data_iso) / (nproc * nsets * nstart * 4))
    phase_space_data_iso = phase_space_data_iso.reshape((nproc, nsets, n_samples_iso, nstart, 4))
