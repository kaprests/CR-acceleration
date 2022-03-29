import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path
from math import pi
from scipy.stats import norm

#DATA_DIR = os.path.dirname(__file__)+'/../Data/'
DATA_DIR = os.path.dirname(__file__)+'/../cluster_dump/randomwalk-data/'
OUT_DIR = '../figs/'
t_max = 20
theta_pi_frac = 0.1
theta = theta_pi_frac*pi
nsets = 100
nstart = 100
nproc = 6
E_inj_exp = 14
z_ax = True
iso_stepsize = False

# Trajectory specific
n_traj_plot = 100
n_steps_plot = 10

D10 = 1.1687289275261758e-05 # E-inj-exp = 10
D14 = 0.11758733574069284 # E-inj-exp = 14
D12 = 0.0011752490186288554 # E-inj-exp = 12

D_dict = {10: D10, 12: D12, 14: D14}
D = D_dict[E_inj_exp]

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

    filename = "randw_"
    if iso_stepsize: filename += "iso-stepsize_"
    if z_ax: filename += "init-z-ax_"
    filename += (
            f"t-max{t_max:.3f}_"
            f"theta-max{theta:.3f}_"
            f"nsets{nsets}_"
            f"nstart{nstart}_"
            f"E-inj-exp{E_inj_exp:.3f}_"
            f"nproc{nproc}"
    )

    filename_iso = (
            f"randw_t-max{t_max:.3f}_"
            f"theta-max{np.pi:.3f}_"
            f"nsets{nsets}_"
            f"nstart{nstart}_"
            f"E-inj-exp{E_inj_exp:.3f}_"
            f"nproc{nproc}"
    )

    # small angle trajectories_
    traj_fname = DATA_DIR + "trajectories_"+filename
    print("============")
    print(traj_fname)
    print("============")
    traj_file = open(traj_fname, 'rb')
    traj_data = np.fromfile(traj_file)
    nsteps = int(len(traj_data)/(nproc*nsets*nstart*4))

    traj_data_proc1 = traj_data[0:int(len(traj_data)/nproc)]
    traj_data_proc1 = traj_data_proc1.reshape((nsets, nsteps, nstart, 4))
    traj_data_proc1_set1 = traj_data_proc1[0]
    traj_data_proc1_set1_particle1 = traj_data_proc1_set1[:, 0, :]
    t, x, y, z = traj_data_proc1_set1_particle1.T

    x = x[0:n_steps_plot]
    y = y[0:n_steps_plot]
    z = z[0:n_steps_plot]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(x, y, z)
    ax.plot(x, y, z, ".")
    plt.title(f"theta_max/pi: {theta_pi_frac}, first {len(x)+1} step positions")
    plt.show()


    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    for i in range(n_traj_plot):
        traj_data_proc1_set1_particlei = traj_data_proc1_set1[:, i, :]
        t, x, y, z = traj_data_proc1_set1_particlei.T
        x = x[0:n_steps_plot]
        y = y[0:n_steps_plot]
        z = z[0:n_steps_plot]
        ax.plot(x, y, z)
        ax.plot(x, y, z, ".")
    plt.title(f"theta_max/pi: {theta_pi_frac}, first {len(x)+1} step positions, for {n_traj_plot} particles")
    plt.show()
