import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

DATA_DIR = './Data/'
# Should instead fix num_steps and sweep opening angle (or cone angle)
num_steps = 10000
basename_iso = 'isotropic_rw'
basename_sa = 'small_angle_rw'

if __name__ == "__main__":
    for num_steps in num_steps_arr:
        # Isotropic rw -- distances
        file = FortranFile(f'{DATA_DIR}{basename_iso}_fdist_{num_steps}', 'r')
        data_iso_dist = file.read_reals()
