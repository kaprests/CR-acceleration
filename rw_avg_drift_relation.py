import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from math import pi
import sys

if __name__ == "__main__":
    ## Isotropic rw -- distances
    print(sys.argv)

    DATA_DIR = './Data/'
    t_max = 101
    theta_pi_fracs = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])
    theta_arr = np.array([pi*f for f in theta_pi_fracs])
    basename_iso = 'isotropic_rw'
    basename_sa = 'small_angle_rw'
    basename_pa = 'pitch_angle_rw'

    iso_dist_avg = np.zeros(len(theta_arr))
    sa_dist_avg = np.zeros(len(theta_arr))
    pa_dist_avg = np.zeros(len(theta_arr))

    for idx, theta in enumerate(theta_arr):
        file = FortranFile(f'{DATA_DIR}{basename_iso}_fdist_tmax{t_max:.3f}', 'r')
        data_iso_dist = file.read_reals()
        iso_dist_avg[idx] = np.average(data_iso_dist)

        file = FortranFile(f'{DATA_DIR}{basename_sa}_fdist_tmax{t_max:.3f}_theta{theta:.3f}')
        data_sa_dist = file.read_reals()
        sa_dist_avg[idx] = np.average(data_sa_dist)

        file = FortranFile(f'{DATA_DIR}{basename_pa}_fdist_tmax{t_max:.3f}_theta{theta:.3f}')
        data_pa_dist = file.read_reals()
        pa_dist_avg[idx] = np.average(data_pa_dist)
    
    rel_avg_drift = sa_dist_avg/iso_dist_avg

#    plt.plot(theta_arr, rel_avg_drift, label="rel drift, sa/iso")
    plt.plot(theta_arr, iso_dist_avg, label="isotropic")
    plt.plot(theta_arr, sa_dist_avg, label="pitch angle, isotropic stepsizes")
    plt.plot(theta_arr, pa_dist_avg, label="pitch angle, adjusted stepsizes")
#    plt.plot(theta_arr, sa_dist_avg*theta_pi_fracs, label="test")
#    plt.plot(theta_arr, [1/(50*theta) for theta in theta_arr], label="K* 1/theta")
    plt.xlabel('theta max')
    plt.ylabel('average drift distance')
    plt.legend()
    plt.show()

    plt.plot(theta_pi_fracs, iso_dist_avg, label="isotropic")
    plt.plot(theta_pi_fracs, sa_dist_avg, label="pitch angle, isotropic stepsizes")
    plt.plot(theta_pi_fracs, pa_dist_avg, label="pitch angle, adjusted stepsizes")
#    plt.plot(theta_pi_fracs, [1/(50*theta) for theta in theta_arr], label="K* 1/theta")
    plt.xlabel('theta_max/pi')
    plt.ylabel('average drift distance')
    plt.legend()
    plt.show()

#    plt.plot(theta_arr, iso_dist_avg, label="isotropic")
#    plt.plot(theta_arr, sa_dist_avg*(pi/100*theta_arr)+0.001*4.5, label="dsa * pi/theta")
#    plt.legend()
#    plt.show()
