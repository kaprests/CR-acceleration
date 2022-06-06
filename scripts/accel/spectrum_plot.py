import numpy as np
import matplotlib.pyplot as plt
from numbers import Number

import os.path
from scipy.optimize import curve_fit

#DATA_DIR = os.path.dirname(__file__) + "/../Data/"
DATA_DIR = os.path.dirname(__file__) + "/../../cluster_data/accel-data/"

#theta_max_pi_frac = 0.0031
theta_max_pi_frac = 1.0
v_shock = 0.01

gamma_shock = 100
n_sets = 100
n_start = 100
E_inj_exp = 10
n_proc = 6
gamma_set = False
z_ax = False


def power_law(x, a, b):
    return a*np.power(x, b)


def spectrum_plot(v_or_gamma, theta_max_pi_frac_array, gamma_set):
    nsets = n_sets
    if gamma_set:
        gamma_shock = v_or_gamma
        basename = f"gamma{gamma_shock:.3f}_"
    else:
        if isinstance(v_or_gamma, Number):
            v_shock = v_or_gamma
            basename = f"vshock{v_shock:.3f}_"
        else:
            # Needs then to be either "injmod1" or "injmod2"
            basename = f"{v_or_gamma}_"
            nsets = 1000
    if z_ax: basename += "init-z-ax_"
    fig, (ax1, ax2) = plt.subplots(1, 2)

    for (i, theta_max_pi_frac) in enumerate(theta_max_pi_frac_array):
        # Load data and plot
        theta = np.pi * theta_max_pi_frac
        N_tot0 = nsets * n_start * n_proc
        filename = (
                f"{basename}"
                f"theta-max{theta:.3f}_"
                f"nsets{nsets}_"
                f"nstart{n_start}_"
                f"E-inj-exp{E_inj_exp:.3f}_"
                f"nproc{n_proc}"
        )
        #spec_fname = DATA_DIR + "spec_prot_"+filename
        #spec_data = np.genfromtxt(spec_fname).T
        #E, m, logE, logm = spec_data
        filename = DATA_DIR + "enumerate_spec_prot_" + filename
        data = np.genfromtxt(filename).T
        E, N0, logE, logE_N0 = data

        N_tot = np.sum(N0)
        N_gte_E = np.zeros(len(E))
        for i, n in enumerate(N0):
            N_gte_E[i] = np.sum(N0[i:])
        #print(N_gte_E[0] / N_tot)
        F = N_gte_E/N_tot

        # E vs F from 'enumerate' data file
        ax1.step(E, N0/(N_tot*E), label=fr"$\theta_{{\mathrm{{max}} }}/\pi = {theta_max_pi_frac}$")

        # E vs F from 'enumerate' data file
        ax2.step(E, E**2 * N0/(N_tot*E), label=fr"$\theta_{{\mathrm{{max}} }}/\pi = {theta_max_pi_frac}$")

    ax1_y0 = (N0/(N_tot*E))[0] * E[0]**2
    ax1.plot(E, ax1_y0/(E**2), label=r"~$1/E^2$")
    ax1.grid(ls="--")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$E$ [fix units]")
    ax1.set_ylabel(r"$F$ [fix units]")
    ax1.legend()

    ax2.grid(ls="--")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$E$ [Fix units]")
    ax2.set_ylabel(r"$F\cdot E^2$ [Fix units]")
    ax2.legend()

    if __name__ == "__main__":
        figname = (
                f"{basename}"
                f"{theta_max_pi_frac_array[0]}_"
                f"nsets{nsets}_"
                f"nstart{n_start}_"
                f"E-inj-exp{E_inj_exp:.3f}_"
                f"nproc{n_proc}"
        )
        plt.savefig(f"spectrum_plot_{figname}.pdf")
        plt.show()



def cross_angle_distribution_plot(v_or_gamma, theta_max_pi_frac_array, gamma_set):
    nsets = n_sets
    if gamma_set:
        gamma_shock = v_or_gamma
        basename = f"gamma{gamma_shock:.3f}_"
    else:
        if isinstance(v_or_gamma, Number):
            v_shock = v_or_gamma
            basename = f"vshock{v_shock:.3f}_"
        else:
            # Needs then to be either "injmod1" or "injmod2"
            basename = f"{v_or_gamma}_"
            nsets = 1000
    if z_ax: basename += "init-z-ax_"
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    for (i, theta_max_pi_frac) in enumerate(theta_max_pi_frac_array):
        # Load data and plot
        theta = np.pi * theta_max_pi_frac
        N_tot0 = nsets * n_start * n_proc
        filename = (
                f"{basename}"
                f"theta-max{theta:.3f}_"
                f"nsets{nsets}_"
                f"nstart{n_start}_"
                f"E-inj-exp{E_inj_exp:.3f}_"
                f"nproc{n_proc}"
        )
        #spec_fname = DATA_DIR + "spec_prot_"+filename
        #spec_data = np.genfromtxt(spec_fname).T
        #E, m, logE, logm = spec_data
        cross_angle_filename = DATA_DIR + "cross-angle-dist_" + filename
        data = np.genfromtxt(cross_angle_filename).T
        theta, initial_down_cross, down_cross, up_cross = data

        if gamma_set:
            aniso_cross_angle_filename = DATA_DIR + "aniso-cross-angle-dist_" + filename
            data = np.genfromtxt(aniso_cross_angle_filename).T
            aniso_theta, aniso_down_cross = data

        ax1.step(theta, initial_down_cross, label=fr"$\theta_{{\mathrm{{max}} }}/\pi = {theta_max_pi_frac}$")
        if gamma_set:
            ax2.set_xlim(right=0.1)
            ax2.step(aniso_theta, aniso_down_cross, label=fr"$\theta_{{\mathrm{{max}} }}/\pi = {theta_max_pi_frac}$")
        else:
            ax2.step(theta, down_cross, label=fr"$\theta_{{\mathrm{{max}} }}/\pi = {theta_max_pi_frac}$")
        ax3.step(theta, up_cross, label=fr"$\theta_{{\mathrm{{max}} }}/\pi = {theta_max_pi_frac}$")
    ax1.set_xlabel(r"$\alpha$ First down crossing")
    ax2.set_xlabel(r"$\alpha$ Ensuing down crossings")
    ax3.set_xlabel(r"$\alpha$ Up crossings")
    ax1.set_ylabel(r"distribution (occurences)")
    ax2.set_ylabel(r"distribution (occurences)")
    ax3.set_ylabel(r"distribution (occurences)")
    ax1.legend(bbox_to_anchor=(1.5, 1.0), loc="lower center", ncol=len(theta_max_pi_frac_array))
    #ax2.legend()
    #ax3.legend()
    if __name__ == "__main__":
        figname = (
                f"{basename}"
                f"nsets{nsets}_"
                f"nstart{n_start}_"
                f"E-inj-exp{E_inj_exp:.3f}_"
                f"nproc{n_proc}"
        )
        plt.savefig(f"cross_angle_dist_plot_{figname}.pdf")
        plt.show()


def aniso_spectrum_cross_angle_plot(v_or_gamma, theta_max_pi_frac_array, gamma_set, aniso=True):
    nsets = n_sets
    if gamma_set:
        gamma_shock = v_or_gamma
        basename = f"gamma{gamma_shock:.3f}_"
    else:
        if isinstance(v_or_gamma, Number):
            v_shock = v_or_gamma
            basename = f"vshock{v_shock:.3f}_"
        else:
            # Needs then to be either "injmod1" or "injmod2"
            basename = f"{v_or_gamma}_"
            nsets = 1000
    if z_ax: basename += "init-z-ax_"
    fig, (ax1, ax2) = plt.subplots(1, 2)

    for (i, theta_max_pi_frac) in enumerate(theta_max_pi_frac_array):
        # Load data and plot
        theta = np.pi * theta_max_pi_frac
        N_tot0 = nsets * n_start * n_proc
        filename = (
                f"{basename}"
                f"theta-max{theta:.3f}_"
                f"nsets{nsets}_"
                f"nstart{n_start}_"
                f"E-inj-exp{E_inj_exp:.3f}_"
                f"nproc{n_proc}"
        )
        #spec_fname = DATA_DIR + "spec_prot_"+filename
        #spec_data = np.genfromtxt(spec_fname).T
        #E, m, logE, logm = spec_data
        spec_filename = DATA_DIR + "enumerate_spec_prot_" + filename
        data = np.genfromtxt(spec_filename).T
        E, N0, logE, logE_N0 = data

        N_tot = np.sum(N0)
        N_gte_E = np.zeros(len(E))
        for i, n in enumerate(N0):
            N_gte_E[i] = np.sum(N0[i:])
        #print(N_gte_E[0] / N_tot)
        F = N_gte_E/N_tot

        if aniso:
            aniso_cross_angle_filename = DATA_DIR + "aniso-cross-angle-dist_" + filename
            data = np.genfromtxt(aniso_cross_angle_filename).T
            aniso_theta, aniso_down_cross = data
            aniso_theta = aniso_theta[2:]
            aniso_down_cross = aniso_down_cross[2:]
        else:
            cross_angle_filename = DATA_DIR + "cross-angle-dist_" + filename
            data = np.genfromtxt(cross_angle_filename).T
            theta, down_cross0, down_cross, up_cross = data

        # E vs F from 'enumerate' data file
        ax1.step(E, N0/(N_tot*E), label=fr"$\theta_{{\mathrm{{max}} }}/\pi = {theta_max_pi_frac}$")

        if aniso:
            ax2.set_xlim(right=0.1)
            ax2.step(aniso_theta, aniso_down_cross/np.sum(aniso_down_cross), label=fr"$\theta_{{\mathrm{{max}} }}/\pi = {theta_max_pi_frac}$")
        else:
            ax2.step(theta, down_cross/np.sum(down_cross), label=fr"$\theta_{{\mathrm{{max}} }}/\pi = {theta_max_pi_frac}$")

    ax1_y0 = (N0/(N_tot*E))[0] * E[0]**2
    ax1.plot(E, ax1_y0/(E**2), label=r"~$1/E^2$")
    ax1.grid(ls="--")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$E$ [fix units]")
    ax1.set_ylabel(r"$F$ [fix units]")
    ax1.legend()

    ax2.grid(ls="--")
    ax2.set_xlabel(r"$\alpha$ Ensuing down crossings")
    ax2.set_ylabel(r"distribution (occurences)")
    ax2.legend()

    if __name__ == "__main__":
        figname = (
                f"{basename}"
                f"{theta_max_pi_frac_array[0]}_"
                f"nsets{nsets}_"
                f"nstart{n_start}_"
                f"E-inj-exp{E_inj_exp:.3f}_"
                f"nproc{n_proc}"
        )
        plt.savefig(f"spectrum_an_cross_angle_plot_{figname}.pdf")
        plt.show()


if __name__ == "__main__":
    ######################
    ### Spectrum plots ###
    ######################

    # Non-relativistic 
    #spectrum_plot(0.01, [1.0, 0.1], False)
    #spectrum_plot("injmod1", [1.0, 0.1], False)
    #spectrum_plot("injmod2", [1.0, 0.1], False)

    ## Ultra-relativistic
    #spectrum_plot(100, [0.0031, 0.006, 0.01], True)
    #spectrum_plot(100, [0.1, 0.5, 1.0], True)

    ######################################
    ### Cross angle distribution plots ###
    ######################################

    # Non-rel
    cross_angle_distribution_plot(0.01, [1.0, 0.1], False)

    # Ultra-rel
    cross_angle_distribution_plot(100, [0.0031, 0.006, 0.01], True)

    # spectrum and angles combines, ultra-rel
    aniso_spectrum_cross_angle_plot(100, [0.0031, 0.006, 0.01], True, aniso=True)
    aniso_spectrum_cross_angle_plot(100, [1.0, 0.5, 0.1], True, aniso=False)
