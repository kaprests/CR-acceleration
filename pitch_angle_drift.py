import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from math import pi
from matplotlib.animation import FuncAnimation
import sys

DATA_DIR = './Data/'
t_max = 101
theta_pi_frac = 1.0
theta = theta_pi_frac*pi
nsets = 10
nstart = 100
stepexp = 2.1

flags = ['--theta-pi-fr', '--tmax', '--stepexp']

print(sys.argv)

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
                                print(f"Invalid argument, using default theta_pi_frac{theta_pi_frac}")
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

    ####################
    ## Isotropic data ##
    ####################
    
    # final distances r(t_max) = sqrt(x**2 + y**2 + z**2)
    #file = FortranFile(DATA_DIR+f'isotropic_rw_fdist_tmax{t_max:.3f}', 'r')
    #file = FortranFile(DATA_DIR+f'pitch_angle_rw_fdist_tmax{t_max:.3f}_theta3.142', 'r')
    #data_iso_dist = file.read_reals()
    
    # final positions x, y, z at t final
    #file = FortranFile(DATA_DIR+f'isotropic_rw_pos_tmax{t_max:.3f}', 'r')
    fdist_iso = FortranFile(DATA_DIR+f'pitch_angle_rw_pos_tmax{t_max:.3f}_theta3.142_stepexp{stepexp:.3f}', 'r')

    # first 10000 events (position and time)
    #file = FortranFile(DATA_DIR+f'isotropic_rw_trajectories_tmax{t_max:.3f}', 'r')
    traj_iso = FortranFile(DATA_DIR+f'pitch_angle_rw_trajectories_tmax{t_max:.3f}_theta3.142_stepexp{stepexp:.3f}', 'r')

    # evenly sampled events (position and time)
    #file = FortranFile(DATA_DIR+f'isotropic_rw_samplepos_tmax{t_max:.3f}', 'r')
    sampos_iso = FortranFile(DATA_DIR+f'pitch_angle_rw_samplepos_tmax{t_max:.3f}_theta3.142_stepexp{stepexp:.3f}', 'r')

    ######################
    ## Small angle data ##
    ######################

    # final distances r(t_max) = sqrt(x**2 + y**2 + z**2)
    #file = FortranFile(DATA_DIR+f'pitch_angle_rw_fdist_tmax{t_max:.3f}_theta{theta:.3f}', 'r')
    #data_sa_dist = file.read_reals()

    # final positions x, y, z at t final
    fdist_pa = FortranFile(DATA_DIR+f'pitch_angle_rw_pos_tmax{t_max:.3f}_theta{theta:.3f}_stepexp{stepexp:.3f}', 'r')

    # first 10000 events (position and time)
    traj_pa = FortranFile(DATA_DIR+f'pitch_angle_rw_trajectories_tmax{t_max:.3f}_theta{theta:.3f}_stepexp{stepexp:.3f}', 'r')

    # evenly sampled events (position and time)
    sampos_pa = FortranFile(DATA_DIR+f'pitch_angle_rw_samplepos_tmax{t_max:.3f}_theta{theta:.3f}_stepexp{stepexp:.3f}', 'r')


    ############################
    ## Plot - Final positions ##
    ############################

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_title("Final positions")

    final_drift_iso = []
    final_drift_pa = []
    for i in range(nsets):
        print(i)
        isopos = fdist_iso.read_reals().reshape(nstart, 3)
        xi = isopos.T[0]
        yi = isopos.T[1]
        zi = isopos.T[2]
        final_drift_iso = np.append(final_drift_iso, xi**2 + yi**2 + zi**2)

        papos = fdist_pa.read_reals().reshape(nstart, 3)
        xp = papos.T[0]
        yp = papos.T[1]
        zp = papos.T[2]
        final_drift_pa = np.append(final_drift_pa, xp**2 + yp**2 + zp**2)
        
        if i == 0:
            ax.scatter(xi, yi, zi, color='green', marker='^', label='isotropic')
            ax.scatter(xp, yp, zp, color='orange', marker='o', label='pitch angle')
        else:
            ax.scatter(xi, yi, zi, color='green', marker='^')
            ax.scatter(xp, yp, zp, color='orange', marker='o')
    avg_drift_iso = np.average(np.sqrt(final_drift_iso))
    avg_drift_pa = np.average(np.sqrt(final_drift_pa))

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend()
    ax.set_title("Final positions")
    plt.show()

    #####################################
    ## Plot - Final drift distribution ##
    #####################################

    plt.hist(final_drift_iso, label='isotropic')
    plt.hist(final_drift_pa, fill=False, label='pitch angle')
    plt.legend()
    plt.title("drift distance distribution")
    plt.show()

    # Print final average distances
    print("Average distances: ")
    print("Average drift distance (isotropic): ", avg_drift_iso)
    print("Average drift distance (pitch angle): ", avg_drift_pa)
    print("Difference: ", avg_drift_pa - avg_drift_iso)
    print()

    ####################################################
    ## Plot - Sampled average drifts vs. sample times ##
    ####################################################

    # Isotropic
    set1_iso = sampos_iso.read_reals()
    num_samples = int(len(set1_iso)/(4*nstart))
    set1_iso = np.reshape(set1_iso, (num_samples, nstart, 4))

    t_sample_iso = set1_iso[:, 1, 3]
    drifts_sampled_iso = np.zeros((num_samples, nsets*nstart))

    for i in range(len(t_sample_iso)):
        x = set1_iso[i, :, 0]
        y = set1_iso[i, :, 1]
        z = set1_iso[i, :, 2]
        n = len(x)
        drifts_sampled_iso[i, 0:len(x)] = x**2 + y**2 + z**2

    for i in range(1, nsets):
        set_iso = sampos_iso.read_reals()
        set_iso = np.reshape(set_iso, (num_samples , nstart, 4))

        for j in range(len(t_sample_iso)):
            x = set_iso[j , :, 0]
            y = set_iso[j , :, 1]
            z = set_iso[j , :, 2]
            n = len(x)
            drifts_sampled_iso[j, i*len(x):(i+1)*len(x)] = x**2 + y**2 + z**2
    avg_drifts_sampled_iso = np.zeros(num_samples)
    for i in range(num_samples):
        avg_drifts_sampled_iso[i] = np.average(np.sqrt(drifts_sampled_iso[i, :]))

    # pitch angle
    set1_pa = sampos_pa.read_reals()
    num_samples = int(len(set1_pa)/(4*nstart))
    set1_pa = np.reshape(set1_pa, (num_samples, nstart, 4))

    t_sample_pa = set1_pa[:, 1, 3]
    drifts_sampled_pa = np.zeros((num_samples, nsets*nstart))

    for i in range(len(t_sample_pa)):
        x = set1_pa[i, :, 0]
        y = set1_pa[i, :, 1]
        z = set1_pa[i, :, 2]
        n = len(x)
        drifts_sampled_pa[i, 0:len(x)] = x**2 + y**2 + z**2

    for i in range(1, nsets):
        set_pa = sampos_pa.read_reals()
        set_pa = np.reshape(set_pa, (num_samples , nstart, 4))

        for j in range(len(t_sample_pa)):
            x = set_pa[j , :, 0]
            y = set_pa[j , :, 1]
            z = set_pa[j , :, 2]
            n = len(x)
            drifts_sampled_pa[j, i*len(x):(i+1)*len(x)] = x**2 + y**2 + z**2
    avg_drifts_sampled_pa = np.zeros(num_samples)
    for i in range(num_samples):
        avg_drifts_sampled_pa[i] = np.average(np.sqrt(drifts_sampled_pa[i, :]))

    plt.plot(t_sample_iso, avg_drifts_sampled_iso, label='isotropic')
    plt.plot(t_sample_pa, avg_drifts_sampled_pa, label='pitch angle')
    plt.title("Average drift distance vs time")
    plt.legend()
    plt.show()

    ###############
    ## animation ##
    ###############

    #fig = plt.figure()
    #ax = fig.add_subplot(projection='3d')
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
    #ax.set_zlabel('z')
    #pos0 = data_sa_trajectories[:, 0, :]
    #x = pos0[0]
    #y = pos0[1]
    #z = pos0[2]
    #scat = ax.scatter(x, y, z, marker='o')

    #def animate(i):
    #    pos = data_sa_trajectories[:, i, :]
    #    x = pos[0]
    #    y = pos[1]
    #    z = pos[2]
    #    scat._offsets3d = (x, y, z)
    #    #ax.scatter(x, y, z, marker='o')

    #anim = FuncAnimation(fig, animate, frames=10000, interval=20, blit=True)
    #plt.show()
