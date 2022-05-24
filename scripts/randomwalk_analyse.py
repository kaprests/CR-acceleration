import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path
import csv
from math import pi
from scipy.stats import norm
from tqdm import trange


###########################
### Settings/parameters ###
###########################
#DATA_DIR = os.path.dirname(__file__)+'/../Data/'
DATA_DIR = os.path.dirname(__file__)+'/../cluster_dump/randomwalk-data/'
OUT_DIR = '../figs/'
t_max = 20
theta_pi_frac = 1.0
theta = theta_pi_frac*pi
nsets = 100
nstart = 100
nproc = 6
E_inj_exp = 10
z_ax = False
iso_stepsize = True


################################################
### Analytical solution diffusion (gaussian) ###
################################################
D10 = 1.1687289275261758e-05 # E-inj-exp = 10
D14 = 0.11758733574069284 # E-inj-exp = 14
D12 = 0.0011752490186288554 # E-inj-exp = 12
D_dict = {10: D10, 12: D12, 14: D14}
D = D_dict[E_inj_exp]
mean = 0
stddev = np.sqrt(D*t_max)
pdf = lambda x: norm.pdf(x, mean, stddev)


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
                                print(f"Using provided theta_pi_frac: {theta_pi_frac}")
                            except:
                                print(f"Invalid argument, using default theta_pi_frac: {theta_pi_frac}    ")
                        elif j == 1:
                            try:
                                t_max = int(val)
                                print(f"Using provided t_max: {t_max}")
                            except:
                                print(f"Invalid argument, using default t_max: {t_max}")
                        elif j == 2:
                            try:
                                stepexp = float(val)
                                print(f"Using provided stepexp: {stepexp}")
                            except:
                                print(f"Invalid argument, using default stepexp: {stepexp}")
            else:
                print(f"Bad flag {flag} provided, no new parameter value set")
    else:
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


def final_positions_data(base_filename):
    """Reads data file containing the final positions and stores in np arrays"""
    fpos_fname = f"{DATA_DIR}fpos_{base_filename}"
    fpos_file = open(fpos_fname, 'rb')
    fpos_data = np.fromfile(fpos_file)
    fpos_data = fpos_data.reshape((nproc*nsets*nstart, 3))
    x_final, y_final, z_final = fpos_data.T
    return x_final, y_final, z_final


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


def final_positions_plot(base_filename, base_filename_iso):
    """Plots the final positions of rw particles

    :base_filename: TODO
    :returns: TODO

    """
    x_final, y_final, z_final = final_positions_data(base_filename)
    x_final_iso, y_final_iso, z_final_iso = final_positions_data(base_filename_iso)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_title("Final positions")
    ax.scatter(x_final, y_final, z_final, marker='o', label=f'theta-max = {theta:.3f}')
    ax.scatter(x_final_iso, y_final_iso, z_final_iso, marker='^', label=f'theta-max = {np.pi:.3f}')
    ax.legend()
    plt.show()


def final_drift_distribution_plot(base_filename, base_filename_iso, n_bins=50):
    """Plots the final drift distribution"""
    x_final, y_final, z_final = final_positions_data(base_filename)
    x_final_iso, y_final_iso, z_final_iso = final_positions_data(base_filename_iso)

    # Total drifted distance
    final_drift_distances = np.sqrt(x_final**2 + y_final**2 + z_final**2)
    final_drift_distances_iso = np.sqrt(x_final_iso**2 + y_final_iso**2 + z_final_iso**2)
    bins = np.linspace(min(final_drift_distances), max(final_drift_distances), n_bins)
    bins_iso = np.linspace(min(final_drift_distances_iso), max(final_drift_distances_iso), n_bins)
    plt.title("Total drift")
    plt.hist(
        final_drift_distances, label='pitch angle', histtype=u'step', color='blue', bins=bins,
        density=True
    )
    plt.hist(
        final_drift_distances_iso, label='isotropic', histtype=u'step', color='red', bins=bins_iso,
        density=True
    )
    plt.show()

    # Drift z-direction
    zbins = np.linspace(min(z_final), max(z_final), n_bins)
    zbins_iso = np.linspace(min(z_final_iso), max(z_final_iso), n_bins)
    z_arr = np.linspace(min(z_final), max(z_final), 1000)
    plt.title("Distribution along z-axis")
    plt.hist(
        z_final, label=f"theta:{theta:.3f}, initial-along-z:{z_ax}", 
        histtype=u'step', color='blue', bins=zbins, density=True
    )
    plt.hist(
        z_final_iso, label=f"theta:{np.pi:.3f}, initial-iso", 
        histtype=u'step', color='red', bins=zbins_iso, density=True
    )
    plt.plot(z_arr, pdf(z_arr), color="green", label="Analytical")    
    plt.legend()    
    plt.show()

    # Drift x-direction
    xbins = np.linspace(min(x_final), max(x_final), n_bins)
    xbins_iso = np.linspace(min(x_final_iso), max(x_final_iso), n_bins)
    x_arr = np.linspace(min(x_final), max(x_final), 1000)
    plt.title("Distribution along x-axis")
    plt.hist(
        x_final, label=f"theta:{theta:.3f}, initial-along-z:{z_ax}", 
        histtype=u'step', color='blue', bins=xbins, density=True
    )
    plt.hist(
        x_final_iso, label=f"theta:{np.pi:.3f}, initial-iso", 
        histtype=u'step', color='red', bins=xbins_iso, density=True
    )
    plt.plot(x_arr, pdf(x_arr), color="green", label="Analytical")    
    plt.legend()
    plt.show()

    # Drift y-direction
    ybins = np.linspace(min(y_final), max(y_final), n_bins)
    ybins_iso = np.linspace(min(y_final_iso), max(y_final_iso), n_bins)
    y_arr = np.linspace(min(y_final), max(y_final), 1000)
    plt.title("Distribution along y-axis")
    plt.hist(
        y_final, label=f"theta:{theta:.3f}, initial-along-z:{z_ax}", 
        histtype=u'step', color='blue', bins=ybins,density=True
    )
    plt.hist(
        y_final_iso, label=f"theta:{np.pi:.3f}, initial-iso", 
        histtype=u'step', color='red', bins=ybins_iso, density=True
    )
    plt.plot(y_arr, pdf(y_arr), color="green", label="Analytical")    
    plt.legend()
    plt.show()


def average_drift_plot(base_filename, base_filename_iso):
    """Plots avg drift vs time

    :base_filename: TODO
    :base_filename_iso: TODO
    :returns: TODO

    """
    t_sampled, avg_drifts_sampled = sample_positions_data(base_filename)
    t_sampled_iso, avg_drifts_sampled_iso = sample_positions_data(base_filename_iso)
    
    plt.plot(t_sampled, avg_drifts_sampled, label=f"theta: {theta}")
    #plt.plot(t_sampled_iso, avg_drifts_sampled_iso, label="isotropic")
    plt.plot(t_sampled_iso, avg_drifts_sampled_iso, ".", label="isotropic")
    plt.legend()
    plt.show()


def D_coeff_estimate(base_filename):
    t_sampled, avg_drifts_sampled = sample_positions_data(base_filename)
    return (avg_drifts_sampled[-1]-avg_drifts_sampled[0])/(t_sampled[-1] - t_sampled[0])/3


def rw_D_coeff_csv():
    print("Writing Dcoeffs to CSV")
    global theta
    global theta_pi_frac
    theta_initial = theta
    theta_pi_frac_intital = theta_pi_frac
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
    theta_pi_frac = theta_pi_frac_intital


if __name__ == "__main__":
    parse_cmd_args()
    filename, filename_iso = construct_base_filename("randw_")

    # Plot functions
    #final_positions_plot(filename, filename_iso)
    #final_drift_distribution_plot(filename, filename_iso)
    #average_drift_plot(filename, filename_iso)
    rw_D_coeff_csv()
