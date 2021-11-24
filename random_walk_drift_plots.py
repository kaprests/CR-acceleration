import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from math import pi
from matplotlib.animation import FuncAnimation

DATA_DIR = './Data/'
t_max = 101
theta_pi_frac = 0.1
theta = theta_pi_frac*pi

if __name__ == "__main__":
    # Isotropic data
    file = FortranFile(DATA_DIR+f'isotropic_rw_fdist_tmax{t_max:.3f}', 'r')
    #file = FortranFile(DATA_DIR+f'pitch_angle_rw_fdist_tmax{t_max:.3f}_theta3.142', 'r')
    data_iso_dist = file.read_reals()
    file = FortranFile(DATA_DIR+f'isotropic_rw_pos_tmax{t_max:.3f}', 'r')
    #file = FortranFile(DATA_DIR+f'pitch_angle_rw_pos_tmax{t_max:.3f}_theta3.142', 'r')
    data_iso_pos = file.read_reals()
    data_iso_pos = np.reshape(data_iso_pos, (int(len(data_iso_pos)/3), 3))
    file = FortranFile(DATA_DIR+f'pitch_angle_rw_trajectories_tmax{t_max:.3f}_theta{3.142}', 'r')
    data_iso_trajectories = file.read_reals()
    data_iso_trajectories = np.reshape(data_iso_trajectories, (4, 10000, 1000))

    # Small angle data
    file = FortranFile(DATA_DIR+f'pitch_angle_rw_fdist_tmax{t_max:.3f}_theta{theta:.3f}', 'r')
    data_sa_dist = file.read_reals()
    file = FortranFile(DATA_DIR+f'pitch_angle_rw_pos_tmax{t_max:.3f}_theta{theta:.3f}', 'r')
    data_sa_pos = file.read_reals()
    data_sa_pos = np.reshape(data_sa_pos, (int(len(data_sa_pos)/3), 3))
    file = FortranFile(DATA_DIR+f'pitch_angle_rw_trajectories_tmax{t_max:.3f}_theta{theta:.3f}', 'r')
    data_sa_trajectories = file.read_reals()
    print(np.size(data_sa_trajectories))
    data_sa_trajectories = np.reshape(data_sa_trajectories, (4, 10000, 1000))

    # Average drift histogram
    plt.hist(data_iso_dist)
    plt.hist(data_sa_dist, fill=False)
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

    print("Average drift distance (isotropic): ", np.average(data_iso_dist))
    print("Average drift distance (isotropic): ", np.average(np.sqrt(xs_iso**2+ys_iso**2+zs_iso**2)))

    print("Average drift distance (pitch angle): ", np.average(np.sqrt(xs_sa**2+ys_sa**2+zs_sa**2)))
    print("Average drift distance (pitch angle): ", np.average(data_sa_dist))

    # intemediate positions
    for i in [0,1,2,3]:
        pos = data_sa_trajectories[:, i, :]

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        xs_sa = pos[0]
        ys_sa = pos[1]
        zs_sa = pos[2]

        ax.scatter(xs_sa, ys_sa, zs_sa, marker='o')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()

    # Average drift over time
    t_array = data_sa_trajectories[3, :, 1]
    avg_drifts_array = np.zeros(len(t_array))

    t_array_iso = data_sa_trajectories[3, :, 1]
    avg_drifts_array_iso = np.zeros(len(t_array))
    for i in range(len(t_array)):
        x = data_sa_trajectories[0, i, :]
        y = data_sa_trajectories[1, i, :]
        z = data_sa_trajectories[2, i, :]
        avg_drifts_array[i] = np.average(np.sqrt(x**2 + y**2 + z**2))

        xi = data_iso_trajectories[0, i, :]
        yi = data_iso_trajectories[1, i, :]
        zi = data_iso_trajectories[2, i, :]
        avg_drifts_array_iso[i] = np.average(np.sqrt(xi**2 + yi**2 + zi**2))

    plt.plot(t_array, avg_drifts_array, label='pitch angle')
    plt.plot(t_array_iso, avg_drifts_array_iso, label='isotropic')
    plt.legend()
    plt.show()

    # animation
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    pos0 = data_sa_trajectories[:, 0, :]
    x = pos0[0]
    y = pos0[1]
    z = pos0[2]
    scat = ax.scatter(x, y, z, marker='o')

    def animate(i):
        pos = data_sa_trajectories[:, i, :]
        x = pos[0]
        y = pos[1]
        z = pos[2]
        scat._offsets3d = (x, y, z)
        #ax.scatter(x, y, z, marker='o')

    anim = FuncAnimation(fig, animate, frames=10000, interval=20, blit=True)
    plt.show()
