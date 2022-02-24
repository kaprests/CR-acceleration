import numpy as np
import matplotlib.pyplot as plt
from math import pi
import sys
from scipy.stats import norm

DATA_DIR = './rwk_data/'
OUT_DIR = './figs/'
t_max = 120
theta_pi_frac = 0.1
theta = theta_pi_frac*pi
nsets = 100
nstart = 10
nproc = 6
E_inj_exp = 10
#stepexp = 2.1
#stepsize = 0
z_ax = True

D = 1.1753625473872295e-05 # For E_inj_exp = 10
mean = 0
stddev = np.sqrt(D*(t_max-100))
pdf = lambda x: norm.pdf(x, mean, stddev)

flags = ['--theta-pi-fr', '--tmax', '--stepexp']

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
                                print(f"Invalid argument, using default theta_pi_frac{theta_pi_frac}    ")
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
    print()

    #####################################
    ### Final particle positions plot ###
    #####################################

    # pitch angle data
    fpos_fname = DATA_DIR
    if z_ax:
        fpos_fname += (
                f"fpos_rwk_init-z-ax_tmax{t_max:.3f}_"
                f"theta-max{theta:.3f}_"
                f"nsets{nsets}_"
                f"nstart{nstart}_"
                f"E_inj_exp{E_inj_exp:.3f}_"
                f"nproc{nproc}"
            )
    else:
        fpos_fname += (
                f"fpos_rwk_tmax{t_max:.3f}_"
                f"theta-max{theta:.3f}_"
                f"nsets{nsets}_"
                f"nstart{nstart}_"
                f"E_inj_exp{E_inj_exp:.3f}_"
                f"nproc{nproc}"
            )
    fpos_file = open(fpos_fname, 'rb')
    fpos_data = np.fromfile(fpos_file)
    
    fpos_data = fpos_data.reshape((nproc*nsets*nstart, 3))
    x_final, y_final, z_final = fpos_data.T

    # Isotropic data (max theta = 2*pi)
    fpos_fname_iso = DATA_DIR
    fpos_fname_iso += (
            f"fpos_rwk_tmax{t_max:.3f}_"
            f"theta-max{np.pi:.3f}_"
            f"nsets{nsets}_"
            f"nstart{nstart}_"
            f"E_inj_exp{E_inj_exp:.3f}_"
            f"nproc{nproc}"
        )
    fpos_file_iso = open(fpos_fname_iso, 'rb')
    fpos_data_iso = np.fromfile(fpos_file_iso)
    
    fpos_data_iso = fpos_data_iso.reshape((nproc*nsets*nstart, 3))
    x_final_iso, y_final_iso, z_final_iso = fpos_data_iso.T

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_title("Final positions")
    ax.scatter(
        x_final, y_final, z_final, marker='o', label=f'theta-max={theta:.3f}, init-z:{z_ax}')
    ax.scatter(x_final_iso, y_final_iso, z_final_iso, marker='^', label=f'theta-max = {np.pi:.3f}')
    ax.legend()
    ax.legend()
    plt.show()

    #############################################
    ###  Final drift distribution (histogram) ###
    #############################################
    
    # Total drift distance
    final_drift_distances = np.sqrt(x_final**2 + y_final**2 + z_final**2)
    final_drift_distances_iso = np.sqrt(x_final_iso**2 + y_final_iso**2 + z_final_iso**2)
    n_bins = 40
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
        histtype=u'step', color='blue', bins=ybins, density=True
    )
    plt.hist(
        y_final_iso, label=f"theta:{np.pi:.3f}, initial-iso", 
        histtype=u'step', color='red', bins=ybins_iso, density=True
    )
    plt.plot(x_arr, pdf(x_arr), color="green", label="Analytical")
    plt.legend()
    plt.show()

    ######################################
    ### Average drift over time - plot ###
    ######################################

    # pitch angle data
    samplepos_fname = DATA_DIR
    if z_ax:
        samplepos_fname += (
                f"samplepos_rwk_init-z-ax_tmax{t_max:.3f}_"
                f"theta-max{theta:.3f}_"
                f"nsets{nsets}_"
                f"nstart{nstart}_"
                f"E_inj_exp{E_inj_exp:.3f}_"
                f"nproc{nproc}"
            )
    else:
        samplepos_fname += (
                f"samplepos_rwk_tmax{t_max:.3f}_"
                f"theta-max{theta:.3f}_"
                f"nsets{nsets}_"
                f"nstart{nstart}_"
                f"E_inj_exp{E_inj_exp:.3f}_"
                f"nproc{nproc}"
            )
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
                x, y, z, t = data.T
                x_samples_seti[proc*nstart:(proc+1)*nstart] = x
                y_samples_seti[proc*nstart:(proc+1)*nstart] = y
                z_samples_seti[proc*nstart:(proc+1)*nstart] = z
                if proc == 0 and seti == 0:
                    t_sampled[sample] = t[0]
            x_samples[seti*len(x_samples_seti):(seti+1)*len(x_samples_seti)] = x_samples_seti
            y_samples[seti*len(y_samples_seti):(seti+1)*len(y_samples_seti)] = y_samples_seti
            z_samples[seti*len(z_samples_seti):(seti+1)*len(z_samples_seti)] = z_samples_seti
        avg_drifts_sampled[sample] = np.average(x_samples**2 + y_samples**2 + z_samples**2)

    # isotropic data
    samplepos_fname_iso = DATA_DIR
    samplepos_fname_iso += (
            f"samplepos_rwk_tmax{t_max:.3f}_"
            f"theta-max{np.pi:.3f}_"
            f"nsets{nsets}_"
            f"nstart{nstart}_"
            f"E_inj_exp{E_inj_exp:.3f}_"
            f"nproc{nproc}"
        )
    samplepos_file_iso = open(samplepos_fname_iso, 'rb')
    samplepos_data_iso = np.fromfile(samplepos_file_iso, dtype='float64')
    n_samples_iso = int(len(samplepos_data_iso) / (nproc * nsets * nstart * 4))
    samplepos_data_iso = samplepos_data_iso.reshape((nproc, nsets, n_samples_iso, nstart, 4))

    avg_drifts_sampled_iso = np.zeros(n_samples)
    t_sampled_iso = np.zeros(n_samples)

    for sample in range(n_samples_iso):
        x_samples = np.zeros(nsets*nstart*nproc)
        y_samples = np.zeros(nsets*nstart*nproc)
        z_samples = np.zeros(nsets*nstart*nproc)
        for seti in range(nsets):
            x_samples_seti = np.zeros(nstart*nproc)
            y_samples_seti = np.zeros(nstart*nproc)
            z_samples_seti = np.zeros(nstart*nproc)
            for proc in range(nproc):
                data = samplepos_data_iso[proc, seti, sample]
                x, y, z, t = data.T
                x_samples_seti[proc*nstart:(proc+1)*nstart] = x
                y_samples_seti[proc*nstart:(proc+1)*nstart] = y
                z_samples_seti[proc*nstart:(proc+1)*nstart] = z
                if proc == 0 and seti == 0:
                    t_sampled_iso[sample] = t[0]
            x_samples[seti*len(x_samples_seti):(seti+1)*len(x_samples_seti)] = x_samples_seti
            y_samples[seti*len(y_samples_seti):(seti+1)*len(y_samples_seti)] = y_samples_seti
            z_samples[seti*len(z_samples_seti):(seti+1)*len(z_samples_seti)] = z_samples_seti
        avg_drifts_sampled_iso[sample] = np.average(x_samples**2 + y_samples**2 + z_samples**2)

    #plt.plot(t_sampled, avg_drifts_sampled, label=f"theta: {theta}")
    #plt.plot(t_sampled_iso, avg_drifts_sampled_iso, label="isotropic")
    #plt.legend()
    #plt.show()

    #print("D': ", (avg_drifts_sampled[-1]-avg_drifts_sampled[0])/(t_sampled[-1] - t_sampled[0])/3)
    #print("D'(iso): ", (avg_drifts_sampled_iso[-1]-avg_drifts_sampled_iso[0])/(t_sampled_iso[-1] - t_sampled_iso[0])/3)
