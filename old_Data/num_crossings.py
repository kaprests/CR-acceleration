import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

DATA_DIR = './Data/'

if __name__ == "__main__":
    # data
    file = FortranFile(DATA_DIR+f'num_crossings_mod_vshock0.010_nsets10_nstart100', 'r')
    file = FortranFile(DATA_DIR+f'num_crossings_mod_gamma10.000_nsets10_nstart100', 'r')
    data = np.sort(file.read_reals())
    N = [0]*len(data)
    n = [0]*len(data)
    for i in range(len(data)):
        N[i] = len([x for x in data if x >= data[i]])
        n[i] = data[i]

    plt.plot(n, N)
    plt.show()
