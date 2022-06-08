import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path
import csv
from math import pi
from scipy.stats import norm
from tqdm import trange


#################
### Constants ###
#################
q_elementary = 0.30282212088 # HEP(?) Loretnz Heaviside units
B0_reg = 1e-6 # regular B -- Gauss
B0_turb = B0_reg/1.0 # turbulent B -- Gauss
m_proton = 938.272e6 # Ev
c = 0.307 # Speed of light, parsec/year


###########################
### Settings/parameters ###
###########################
#DATA_DIR = os.path.dirname(__file__)+'/../Data/'
DATA_DIR = os.path.dirname(__file__)+'/../../cluster_data/randomwalk-data/'
OUT_DIR = f"{os.path.dirname(__file__)}/figs/"
t_max = 20
theta_pi_frac = 1.0
theta = theta_pi_frac*pi
nsets = 100
nstart = 100
nproc = 6
E_inj_exp = 10
z_ax = False
iso_stepsize = True
save_figs = False


################################################
### Analytical solution diffusion (gaussian) ###
################################################
def v_particle(E, m):
    """
    Comute the velocity of a particle from given energy

    Units where c = 1
    """
    return np.sqrt(E**2 - m**2)/E


def R_larmor_year(E, B=B0_turb):
    """
    Larmor radius from energy in simulation units ([t] = [l] = year, [v] = 1)
    
    Lightspeed c = 1. 
    """
    return 3.523e-21*E/B


def D_coeff_year(E):
    """
    Computes the diffusion coefficient corresponding to the mean free path corresponding 
    to the given particle energy, 
    for a random walk with isotropic scattering. 
    
    Simulation units
    """
    return v_particle(E, m_proton)*R_larmor_year(E)/3


def R_larmor_parsec(E, B=B0_turb):
    """
    Larmor radius in parsec. 
    """
    R_year = R_larmor_year(E, B=B)
    R_parsec = R_year * c
    return R_parsec


def D_coeff_pc2_yr(E):
    """
    Diffusion coefficient in units pc^2/year. 
    """
    v = v_particle(E, m_proton) * c # particle velocity in parsec per year
    return v * R_larmor_parsec(E)/3



#D10 = 1.1687289275261758e-05 # E-inj-exp = 10
#D14 = 0.11758733574069284 # E-inj-exp = 14
#D12 = 0.0011752490186288554 # E-inj-exp = 12
D10 = D_coeff_year(1e10)
D12 = D_coeff_year(1e12)
D14 = D_coeff_year(1e14)
D_dict = {10: D10, 12: D12, 14: D14}
D = D_dict[E_inj_exp]
mean = 0
stddev = np.sqrt(D*t_max)
pdf = lambda x: norm.pdf(x, mean, stddev)


#######################
### setup functions ###
#######################


def parse_cmd_args():
    """Parses cmd line arguments
    :returns: _
    """
    global theta_pi_frac
    global theta
    global t_max
    global stepexp
    flags = ['--theta-pi-fr', '--tmax', '--stepexp']
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
                                if __name__ == "__main__":
                                    print(f"Using provided theta_pi_frac: {theta_pi_frac}")
                            except:
                                if __name__ == "__main__":
                                    print(f"Invalid argument, using default theta_pi_frac: {theta_pi_frac}    ")
                        elif j == 1:
                            try:
                                t_max = int(val)
                                if __name__ == "__main__":
                                    print(f"Using provided t_max: {t_max}")
                            except:
                                if __name__ == "__main__":
                                    print(f"Invalid argument, using default t_max: {t_max}")
                        elif j == 2:
                            try:
                                stepexp = float(val)
                                if __name__ == "__main__":
                                    print(f"Using provided stepexp: {stepexp}")
                            except:
                                if __name__ == "__main__":
                                    print(f"Invalid argument, using default stepexp: {stepexp}")
            else:
                if __name__ == "__main__":
                    print(f"Bad flag {flag} provided, no new parameter value set")
    else:
        if __name__ == "__main__":
            print("Odd number of flags+parameters, using default values")


def construct_base_filename(base):
    """Constructs filename

    :base: TODO
    :returns: TODO

    """
    filename = base
    filename_iso = base
    if iso_stepsize: filename += "iso-stepsize_"
    if z_ax: filename += "init-z-ax_"
    filename += (
            f"t-max{t_max:.3f}_"
            f"theta-max{theta:.3f}_"
            f"nsets{nsets}_"
            f"nstart{nstart}_"
            f"E-inj-exp{E_inj_exp:.3f}_"
            f"nproc{nproc}"
    )

    filename_iso += (
            f"iso-stepsize_t-max{t_max:.3f}_"
            f"theta-max{np.pi:.3f}_"
            f"nsets{nsets}_"
            f"nstart{nstart}_"
            f"E-inj-exp{E_inj_exp:.3f}_"
            f"nproc{nproc}"
    )
    return filename, filename_iso


###################################
### Data reading and processing ###
###################################


def final_positions_data(base_filename):
    """Reads data file containing the final positions and stores in np arrays"""
    fpos_fname = f"{DATA_DIR}fpos_{base_filename}"
    fpos_file = open(fpos_fname, 'rb')
    fpos_data = np.fromfile(fpos_file)
    fpos_data = fpos_data.reshape((nproc*nsets*nstart, 3))
    x_final, y_final, z_final = fpos_data.T
    return x_final, y_final, z_final


# Rewrite: try using some reshaping instead of the foor loops
def sample_positions_data(base_filename):
    """Reads data file containing the sampled positions and stores in np arrays"""
    samplepos_fname = f"{DATA_DIR}samplepos_{base_filename}"
    samplepos_file = open(samplepos_fname, 'rb')
    samplepos_data = np.fromfile(samplepos_file, dtype='float64')
    n_samples = int(len(samplepos_data) / (nproc * nsets * nstart * 4))
    samplepos_data = samplepos_data.reshape((nproc, nsets, n_samples, nstart, 4))

    avg_drifts_sampled = np.zeros(n_samples)
    t_sampled = np.zeros(n_samples)

    for sample in range(n_samples):
        x_samples = np.zeros(nsets*nstart*nproc)
        y_samples = np.zeros(nsets*nstart*nproc)
        z_samples = np.zeros(nsets*nstart*nproc)
        for seti in range(nsets):
            x_samples_seti = np.zeros(nstart*nproc)
            y_samples_seti = np.zeros(nstart*nproc)
            z_samples_seti = np.zeros(nstart*nproc)
            for proc in range(nproc):
                data = samplepos_data[proc, seti, sample]
                t, x, y, z = data.T
                x_samples_seti[proc*nstart:(proc+1)*nstart] = x
                y_samples_seti[proc*nstart:(proc+1)*nstart] = y
                z_samples_seti[proc*nstart:(proc+1)*nstart] = z
                if proc == 0 and seti == 0:
                    t_sampled[sample] = t[0]
            x_samples[seti*len(x_samples_seti):(seti+1)*len(x_samples_seti)] = x_samples_seti
            y_samples[seti*len(y_samples_seti):(seti+1)*len(y_samples_seti)] = y_samples_seti
            z_samples[seti*len(z_samples_seti):(seti+1)*len(z_samples_seti)] = z_samples_seti
        avg_drifts_sampled[sample] = np.average(x_samples**2 + y_samples**2 + z_samples**2)

    t_sampled = t_sampled[:-int(n_samples/10)]
    avg_drifts_sampled = avg_drifts_sampled[:-int(n_samples/10)]
    return t_sampled, avg_drifts_sampled


#############
### Plots ###
#############


#Non production
def final_positions_plot(theta_pi_frac_arr, iso=None):
    """Plots the final positions of rw particles

    :base_filename: TODO
    :returns: TODO

    """
    global theta
    global iso_stepsize
    iso_stepsize_initial = iso_stepsize
    theta_initial = theta
    if iso != None:
        iso_stepsize = iso
    for (i, theta_pi_frac) in enumerate(theta_pi_frac_arr):
        theta = theta_pi_frac * np.pi
        filename, filename_iso = construct_base_filename("randw_")
        x_final, y_final, z_final = final_positions_data(filename)
        x_final_iso, y_final_iso, z_final_iso = final_positions_data(filename_iso)

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_title(fr"Final positions")
        ax.scatter(x_final, y_final, z_final, marker='o', label=fr'$\theta_{{\mathrm{{max}}}}/\pi={theta_pi_frac}$')
        if i == 0:
            ax.scatter(x_final_iso, y_final_iso, z_final_iso, marker='^', label=fr'$\theta_{{\mathrm{{max}}}}/\pi =1.0$')
    theta = theta_initial
    iso_stepsize = iso_stepsize_initial
    ax.legend()
    if __name__ == "__main__":
        plt.show()


def total_final_drift_distribution_plot(theta_pi_frac_arr, n_bins=50, iso=None, datafname=None):
    """Plot of total drift distribution from origin"""
    global iso_stepsize
    iso_stepsize_initial = iso_stepsize
    if iso != None:
        iso_stepsize = iso    
    colors = ["black", "green", "blue", "orange", "yellow"]
    global theta
    theta_initial = theta

    #_, filename_iso = construct_base_filename("randw_")
    #x_final_iso, y_final_iso, z_final_iso = final_positions_data(filename_iso)
    #final_drift_distances_iso = np.sqrt(x_final_iso**2 + y_final_iso**2 + z_final_iso**2)
    #bins_iso = np.linspace(min(final_drift_distances_iso), max(final_drift_distances_iso), n_bins)
    #plt.hist(
    #    final_drift_distances_iso, histtype=u'step', color='red', bins=bins_iso,
    #    density=True, label=fr'$\theta_{{\mathrm{{max}}}}/\pi=1.0$'
    #)

    theta_pi_frac_arr = sorted(theta_pi_frac_arr, reverse=True)
    for (theta_pi_frac, color) in zip(theta_pi_frac_arr, colors): # not global theta!
        theta = theta_pi_frac * np.pi
        filename, _ = construct_base_filename("randw_")
        x_final, y_final, z_final = final_positions_data(filename)

        # Total drifted distance
        final_drift_distances = np.sqrt(x_final**2 + y_final**2 + z_final**2)
        bins = np.linspace(min(final_drift_distances), max(final_drift_distances), n_bins)
        plt.hist(
            final_drift_distances, histtype=u'step', color=color, bins=bins,
            density=True, label=fr'$\theta_{{\mathrm{{max}}}}/\pi={theta_pi_frac}$'
        )
    theta = theta_initial
    iso_stepsize = iso_stepsize_initial

    #plt.title(fr"Total drift")
    plt.xlabel("x [fix units]")
    plt.ylabel("y [fix units]")
    plt.legend()
    pdf_name = "final_total_drift_selected_theta"
    if iso_stepsize: pdf_name += "_no-step-corr"
    if z_ax: pdf_name += "_init-z"
    if __name__ == "__main__":
        plt.show()


def final_drift_distribution_plot(theta_pi_frac_arr, n_bins=50, legend=True, iso=None):
    """Plots the final drift distribution along each axis"""
    # Get data
    colors = ["black", "green", "blue", "orange", "yellow"]
    global theta
    global iso_stepsize
    theta_initial = theta
    iso_stepsize_initial = iso_stepsize
    if iso != None:
        iso_stepsize = iso
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
    #fig.suptitle("Distribution along each spatial axis")
    axes = [ax1, ax2, ax3]
    for (i, theta_pi_frac) in enumerate(theta_pi_frac_arr):
        theta = theta_pi_frac * np.pi
        filename, filename_iso = construct_base_filename("randw_")
        x_final, y_final, z_final = final_positions_data(filename)
        x_final_iso, y_final_iso, z_final_iso = final_positions_data(filename_iso)

        # Drift x-direction
        xbins = np.linspace(min(x_final), max(x_final), n_bins)
        xbins_iso = np.linspace(min(x_final_iso), max(x_final_iso), n_bins)
        x_arr = np.linspace(min(x_final), max(x_final), 1000)
        # Drift y-direction
        ybins = np.linspace(min(y_final), max(y_final), n_bins)
        ybins_iso = np.linspace(min(y_final_iso), max(y_final_iso), n_bins)
        y_arr = np.linspace(min(y_final), max(y_final), 1000)
        zbins = np.linspace(min(z_final), max(z_final), n_bins)
        zbins_iso = np.linspace(min(z_final_iso), max(z_final_iso), n_bins)
        z_arr = np.linspace(min(z_final), max(z_final), 1000)

        axis_labels = ['x', 'y', 'z']
        ax_bins_array = [xbins, ybins, zbins]
        ax_bins_array_iso = [xbins_iso, ybins_iso, zbins_iso]
        ax_arrays = [x_arr, y_arr, z_arr]
        ax_final_arrays = [x_final, y_final, z_final]
        ax_final_arrays_iso = [x_final_iso, y_final_iso, z_final_iso]

        for (j, axis_label) in enumerate(axis_labels):
            ax_array = ax_arrays[j]
            ax_final_array = ax_final_arrays[j]
            if i == 0:
                axes[j].plot(ax_array, pdf(ax_array), linestyle=(0, (5, 1)), color="gray", label="Target distribution")
            hist = axes[j].hist(
                ax_final_array, 
                label=fr'$\theta_{{\mathrm{{max}}}}/\pi={theta_pi_frac}$', 
                histtype=u'step', color=colors[i], bins=ax_bins_array[j], density=True
            )
            if i == 0:
                axes[j].set_xlabel(f"{axis_labels[j]}-axis [fix units]")
                axes[j].set_ylim(top=max(hist[0])*1.05)
                if j == 0:
                    axes[j].set_ylabel(f"Distribution [fix units]")
    theta = theta_initial
    iso_stepsize = iso_stepsize_initial
    handles, labels = ax1.get_legend_handles_labels()
    if legend:
        if __name__ == "__main__":
            fig.legend(handles, labels, loc='upper center', ncol=len(theta_pi_frac_arr)+1, mode="expand")
        else: 
            # tikzplotlib friendly legend
            ax2.legend(bbox_to_anchor=(0.5, 1.1), loc="lower center", ncol=len(theta_pi_frac_arr)+1)
    if __name__ == "__main__":
        plt.show()


def average_drift_plot(theta_pi_frac_arr, ax=None):
    """Plots avg drift vs time

    :base_filename: TODO
    :base_filename_iso: TODO
    :returns: TODO

    """
    global theta
    theta_initial = theta
    target_plotted = False
    for (i, theta_pi_frac) in enumerate(theta_pi_frac_arr):
        theta = theta_pi_frac * np.pi
        filename, filename_iso = construct_base_filename("randw_")
        t_sampled, avg_drifts_sampled = sample_positions_data(filename)

        if theta_pi_frac == 1.0:
            continue
        if not target_plotted:
            t_sampled_iso, avg_drifts_sampled_iso = sample_positions_data(filename_iso)
            if ax:
                ax.plot(t_sampled_iso, avg_drifts_sampled_iso, color="gray", label="Target drift") # replace with theoretical
                ax.plot(t_sampled_iso[0::20], avg_drifts_sampled_iso[0::20], "+", 
                        color="red", markersize=5,label=fr"$\theta_{{\mathrm{{max}}}}/\pi=1.0$")
            else:
                plt.plot(t_sampled_iso, avg_drifts_sampled_iso, color="gray", label="Target drift") # replace with theoretical
                plt.plot(t_sampled_iso[0::20], avg_drifts_sampled_iso[0::20], "+", 
                        color="red", markersize=5,label=fr"$\theta_{{\mathrm{{max}}}}/\pi=1.0$")
            target_plotted = True
        if ax:
            ax.plot(t_sampled[0::10], avg_drifts_sampled[0::10], linestyle="-", label=fr"$\theta_{{\mathrm{{max}} }}/\pi={theta_pi_frac}$")
        else:
            plt.plot(t_sampled[0::10], avg_drifts_sampled[0::10], linestyle="-", label=fr"$\theta_{{\mathrm{{max}} }}/\pi={theta_pi_frac}$")

        #plt.plot(t_sampled_iso[0::10], avg_drifts_sampled_iso[0::10], "+", color="red", label=fr"$\theta_{{\mathrm{{max}} }}/\pi=1.0$")
    if ax:
        ax.set_xlabel("x [add untis]")
        ax.set_ylabel("y [add untis]")
        ax.legend()
    else:
        plt.xlabel("x [add untis]")
        plt.ylabel("y [add untis]")
        plt.legend()
    if __name__ == "__main__" and not ax:
        plt.show()


########################
### Computiing stuff ###
########################


def D_coeff_estimate(base_filename):
    t_sampled, avg_drifts_sampled = sample_positions_data(base_filename)
    return (avg_drifts_sampled[-1]-avg_drifts_sampled[0])/(t_sampled[-1] - t_sampled[0])/3


def rw_D_coeff_csv():
    print("Writing Dcoeffs to CSV")
    global theta_pi_frac
    theta_initial = theta
    theta_pi_frac_initial = theta_pi_frac
    with open('rw_dcoeff.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow([f"{0.0:.3f}", 'inf'])
        theta0_arr = [0.01, 0.1]
        for theta0 in theta0_arr:
            for i in trange(9):
                print(f"{theta0+i*theta0:.3f}")
                theta_pi_frac = float(f"{theta0+i*theta0:.3f}")
                theta = theta_pi_frac * np.pi
                filename, filename_iso = construct_base_filename("randw_")
                writer.writerow([f"{theta_pi_frac:.3f}", D_coeff_estimate(filename)])
        writer.writerow([f"{1.0:.3f}", D_coeff_estimate(filename_iso)])
    theta = theta_initial
    theta_pi_frac = theta_pi_frac_initial


if __name__ == "__main__":
    # Setup
    parse_cmd_args() # Parse optional cmd args -- overrides default parameter values
    #theta_pi_frac_selected = [0.05, 0.1, 0.5]
    theta_pi_frac_selected = [0.1, 0.5, 1.0]
    theta_pi_frac_all = [0.05, 0.1, 0.5]
    filename, filename_iso = construct_base_filename("randw_") # construct base filename of the correct data file(s)

    # Final positions scatter 3D plot
    #final_positions_plot([theta_pi_frac])

    # Total drift
    #total_final_drift_distribution_plot([theta_pi_frac]) # Single theta 
    #total_final_drift_distribution_plot(theta_pi_frac_selected) # Multiple theta

    # Drift along axis
    #final_drift_distribution_plot(theta_pi_frac_selected)
    #final_drift_distribution_plot([theta_pi_frac])

    # average drift over time
    #average_drift_plot(filename, filename_iso)

    # Compute corresponding diffusion coefficients
    #rw_D_coeff_csv()

    #######################
    # Test for production #
    #######################

    #total_final_drift_distribution_plot([1.0, 0.7, 0.5, 0.3])
    #total_final_drift_distribution_plot([0.1, 0.07, 0.05, 0.03])

    #final_drift_distribution_plot([1.0, 0.5])
    #final_drift_distribution_plot([0.05, 0.01], legend=False)

    #fig, (ax1, ax2) = plt.subplots(1, 2)
    #average_drift_plot([0.7, 0.5, 0.3], ax=ax1)
    #average_drift_plot([0.3, 0.2, 0.1], ax=ax2)
    #plt.show()
    #average_drift_plot([0.3], ax=ax1)
    #average_drift_plot([0.1], ax=ax2)
    #plt.show()

    #total_final_drift_distribution_plot([1.0, 0.5, 0.1], iso=False)
    #plt.show()
    #final_drift_distribution_plot([0.5], iso=False)
    #final_drift_distribution_plot([0.1], iso=False)
    #plt.show()

    print(D_coeff_year(1e10)/D10)
    print(D_coeff_year(1e12)/D12)
    print(D_coeff_year(1e14)/D14)

    print(D_coeff_pc2_yr(1e10))
    print(D_coeff_year(1e10))
