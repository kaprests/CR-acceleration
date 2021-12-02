import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from math import pi
from matplotlib.animation import FuncAnimation
import sys

DATA_DIR = './Data/'
OUT_DIR = './figs/'
nsets = 100
nstart = 100
stepexp = 2.1
injmod = 2
rw_model = 'iso'
#rw_model = 'pas'

#flags = ['--theta-pi-fr', '--tmax', '--stepexp']
#
#print(sys.argv)
#
if __name__ == "__main__":
#    if (len(sys.argv)-1)%2 == 0:
#        args = sys.argv[1:]
#        for i in range(0, len(args), 2):
#            flag = args[i]
#            val = args[i+1]
#
#            if flag in flags:
#                for j, f in enumerate(flags):
#                    if flag == f:
#                        if j == 0:
#                            try:
#                                theta_pi_frac = float(val)
#                                theta = theta_pi_frac * pi
#                                print(f"Using provided theta_pi_frac{theta_pi_frac}")
#                            except:
#                                print(f"Invalid argument, using default theta_pi_frac{theta_pi_frac}")
#                        elif j == 1:
#                            try:
#                                t_max = int(val)
#                                print(f"Using provided t_max{t_max}")
#                            except:
#                                print(f"Invalid argument, using default t_max{t_max}")
#                        elif j == 2:
#                            try:
#                                stepexp = float(val)
#                                print(f"Using provided stepexp{stepexp}")
#                            except:
#                                print(f"Invalid argument, using default stepexp{stepexp}")
#            else:
#                print(f"Bad flag {flag} provided, no new parameter value set")
#    else:
#        print("Odd number of flags+parameters, using default values")

    ##########
    ## data ##
    ##########

    # first 10000 events (position and time)
    traj = FortranFile(DATA_DIR+f'trajectories_stepexp{stepexp:.3f}_{rw_model}_injmod{injmod}_nsets{nsets}_nstart{nstart}', 'r')
    
    ###############
    ## animation ##
    ###############

    #fig = plt.figure()
    #ax = fig.add_subplot(projection='3d')
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
    #ax.set_zlabel('z')

    num_steps = 10000
    x_arr = np.zeros((num_steps, nsets*nstart))
    y_arr = np.zeros((num_steps, nsets*nstart))
    z_arr = np.zeros((num_steps, nsets*nstart))

    for i in range(nsets):
        dset = traj.read_reals()
        dset = np.reshape(dset, (num_steps, nstart, 4))

        for j in range(num_steps):
            xj = dset[j, :, 0]
            yj = dset[j, :, 1]
            zj = dset[j, :, 2]
            n = len(xj)
            x_arr[j, i*n:(i+1)*n] = xj
            y_arr[j, i*n:(i+1)*n] = yj
            z_arr[j, i*n:(i+1)*n] = zj

    for i in range(30):
        x = x_arr[i, :]
        y = y_arr[i, :]
        z = z_arr[i, :]
        #scat = ax.scatter(x, y, z, marker='o')
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(x, y, z)
        ax.set_title(f'{i+1}/30 -- rw_model: {rw_model}')
        plt.savefig(f"{OUT_DIR}{rw_model}_with_shock_positions{i}.pdf")
        plt.show()

        x_arr = x_arr[1:, :]
        y_arr = y_arr[1:, :]
        z_arr = z_arr[1:, :]


    def animate(i):
        # Called every time step
        xi = x_arr[i, :]
        yi = y_arr[i, :]
        zi = z_arr[i, :]
        scat._offsets3d = (xi, yi, zi)


    #anim = FuncAnimation(fig, animate, frames=10000, interval=20, blit=True)
    #plt.show()
