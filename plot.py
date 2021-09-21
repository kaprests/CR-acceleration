import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

DATA_DIR = './Data/'
RAW_DATA_DIR = './RawData/'


def datapath(n_sets, n_start, version, v_shock=-1, dir=DATA_DIR, base='spec_prot'):
    if v_shock >0:
        vsh = '{:.3f}'.format(v_shock)
        path = f"{dir}{base}_{version}_vshock{vsh}_nsets{n_sets}_nstart{n_start}"
    else:
        path = f"{dir}{base}_{version}_nsets{n_sets}_nstart{n_start}"
    return path


def bin_data(data, nbins):
    """Bins data in array data

    :data: array with data to be binned
    :nbins: Number of bins
    :returns: TODO

    """
    # Everything log10
    data = np.log10(data)
    lower_lim = min(data)
    upper_lim = max(data)
    bins = np.linspace(lower_lim, upper_lim, nbins)
    digitized = np.digitize(data, bins) # Bin number of each element in data
    binsize = bins[1] - bins[0]
    bin_halfsize = binsize/2
    bin_middle = [bin_halfsize + bins[i] for i in range(len(bins)-1)]
    bin_middle = np.array(bin_middle)
    bin_counts = [np.count_nonzero(digitized == i) for i in range(1, nbins)]
    bin_counts = np.array(bin_counts)
    N_gte_En = np.array([np.sum(bin_counts[i:]) for i in range(len(bin_counts))])
    
    return bins, np.array(bin_middle), N_gte_En

if __name__ == "__main__":
    ### config ###
    n_start = 100
    n_sets = 10
    version = "mod"
    nbins = 10 # for binning raw data
    v_shock = 0.100

    ### pre-binned data
    path = datapath(n_sets, n_start, version)
    data = np.fromfile(path, sep=' ')
    data = data.reshape(int(len(data)/4), 4).T
    E, m, logE, logm = data
    plt.step(logE, logm)
    plt.show()
    
    ### raw data ###
    raw_path = datapath(10, 100, 'mod', dir=RAW_DATA_DIR, base='exit_energies', v_shock=0.100)
    raw_path = raw_path + '.dat'
    file = FortranFile(raw_path, 'r')
    data = file.read_reals()
    log_bins, log10_bin_middle, N_gte_En = bin_data(data, nbins)

    bin_energies = 10**(log10_bin_middle)
    m = N_gte_En * bin_energies**2 # should be ~ constant for v_shock << c WRONG m somehow
    
    plt.step(np.log10(bin_energies), np.log10(m))
    plt.show()
