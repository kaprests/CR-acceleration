import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from math import pi


DATA_DIR = './Data/'
OUT_DIR = './figs/'
nsets = 100
nstart = 100
injmod = 2
stepexp = 2.1
rw_model = 'iso'
#rw_model = 'pas'
n_cross = 10


if __name__ == "__main__":
    cross1 = FortranFile(DATA_DIR+f'cross_angles_stepexp{stepexp:.3f}_{rw_model}_injmod{injmod}_nsets{nsets}_nstart{nstart}', 'r').read_reals()
    print(np.sum(cross1))


    cross1 = FortranFile(DATA_DIR+f'cross_angles_stepexp{stepexp:.3f}_{rw_model}_injmod{injmod}_nsets{nsets}_nstart{nstart}', 'r').read_reals()
    print(np.sum(cross1))
