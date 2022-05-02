import numpy as np
import matplotlib.pyplot as plt
import os.path

DATA_DIR = os.path.dirname(__file__) + "/../Data/"

theta_max_pi_frac = 1.0
v_shock = 0.01
gamma_shock = 100
n_sets = 1000
n_start = 10
E_inj_exp = 10
n_proc = 2
gamma_set = False
z_ax = False


if __name__ == "__main__":
    theta = np.pi * theta_max_pi_frac
    N_tot = n_sets * n_sets * n_proc

    if gamma_set:
        filename = f"gamma{gamma_shock:.3f}_"
    else:
        filename = f"vshock{v_shock:.3f}_"
    if z_ax: filename += "init-z-ax_"
    filename += (
            f"theta-max{theta:.3f}_"
            f"nsets{n_sets}_"
            f"nstart{n_start}_"
            f"E-inj-exp{E_inj_exp:.3f}_"
            f"nproc{n_proc}"
    )

    # small angle trajectories_
    spec_fname = DATA_DIR + "spec_prot_"+filename
    spec_data = np.genfromtxt(spec_fname).T
    E, m, logE, logm = spec_data
    
    enumerate_fname = DATA_DIR + "enumerate_spec_prot_" + filename
    enum_data = np.genfromtxt(enumerate_fname).T
    E_enum, N0, logE_enum, logE_N0 = enum_data

    if E != E_enum:
        print("OBS: E != E_enum")

    N_gte_E = np.zeros(len(E))
    for i, n in enumerate(N0):
        N_gte_E[i] = np.sum(N0[i:])

    #assert(N_gte_E[-1] == N0[-1])
    print(len(N_gte_E))
    print(len(N0))
    print(len(E))
    print(N_gte_E[-1])
    print(N0[-1])

    assert(N_gte_E[0]/N_tot == 1)
    F_enum = N_gte_E/N_tot
    
    # Original
    plt.title("Original plot")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("E")
    plt.ylabel("E**2 * F")
    plt.step(E, m)
    plt.show()
    
    # E vs F (F = m/E**2)
    plt.title("E vs F")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("E")
    plt.ylabel("F")
    plt.step(E, m/E**2)
    plt.show()

    # E vs F from 'enumerate' data file
    plt.title("E vs F (new data)")
    plt.xscale("log")
    plt.xscale("log")
    plt.xlabel("E")
    plt.ylabel("F")
    plt.step(E_enum, F_enum)
    plt.show()
