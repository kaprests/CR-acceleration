import numpy as np
import matplotlib.pyplot as plt
from math import pi
import sys
import os.path

DATA_DIR = os.path.dirname(__file__)+'/../Data/'
OUT_DIR = './figs/'
theta_pi_frac = 1.0
theta = theta_pi_frac*pi
nsets = 500
nstart = 100
nproc = 6
injmod = 1
num_phase_log = 5
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
    print(phase_space_data.size)
    print(nsets*nstart*2*num_phase_log)
    #data_array_size = nstart*2*n_log
    #set_size = data_array_size * nsets
    #total_data_size = set_size * nproc
