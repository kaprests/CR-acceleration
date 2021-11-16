import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

DATA_DIR = './Data/'
num_steps = 1000

if __name__ == "__main__":
    # Isotropic data
    file = FortranFile(DATA_DIR+f'isotropic_rw_fdist_{num_steps}', 'r')
    data_iso_dist = file.read_reals()
    file = FortranFile(DATA_DIR+f'isotropic_rw_pos_{num_steps}', 'r')
    data_iso_pos = file.read_reals()
    data_iso_pos = np.reshape(data_iso_pos, (int(len(data_iso_pos)/3), 3))

    # Small angle data
    file = FortranFile(DATA_DIR+f'small_angle_rw_fdist_nsteps{num_steps}_theta3.142', 'r')
    data_sa_dist = file.read_reals()
    file = FortranFile(DATA_DIR+f'small_angle_rw_pos_nsteps{num_steps}_theta3.142', 'r')
    data_sa_pos = file.read_reals()
    data_sa_pos = np.reshape(data_sa_pos, (int(len(data_sa_pos)/3), 3))

    # Average drift histogram
    plt.hist(data_iso_dist)
    plt.hist(data_sa_dist)
    plt.show()

    # Final positions - isotropic_rw_pos
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    xs_iso = data_iso_pos.T[0]
    ys_iso = data_iso_pos.T[1]
    zs_iso = data_iso_pos.T[2]

    xs_sa = data_sa_pos.T[0]
    ys_sa = data_sa_pos.T[1]
    zs_sa = data_sa_pos.T[2]

    ax.scatter(xs_iso, ys_iso, zs_iso, marker='o')
    ax.scatter(xs_sa, ys_sa, zs_sa, marker='o')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

    print("Average drift distance (iso): ", np.average(data_iso_dist))
    print("Average drift distance (sa): ", np.average(data_sa_dist))
