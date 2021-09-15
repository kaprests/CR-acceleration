import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = './Data/'


def datapath(n_sets, n_start, version, v_shock=-1, dir=DATA_DIR):
    if v_shock >0:
        path = f"{dir}spec_prot_{version}_vshock{v_shock}_nsets{n_sets}_nstart{n_start}"
    else:
        path = f"{dir}spec_prot_{version}_nsets{n_sets}_nstart{n_start}"
    return path


if __name__ == "__main__":
    ### config ###
    n_start = 100
    n_sets = 10
    v_shock = -1
    version = "mod"

    path = datapath(n_sets, n_start, version)
    data = np.fromfile(path, sep=' ')
    data = data.reshape(int(len(data)/4), 4).T
    E, m, logE, logm = data
    
    plt.step(logE, logm)
    plt.show()
