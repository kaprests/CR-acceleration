import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import os.path


# Input file
in_file = f"{os.path.dirname(__file__)}/rw_dcoeff.csv"
# Output directory
OUT_DIR = f"{os.path.dirname(__file__)}/../../../cr-acceleration-report/data/"


def save_table(columns, header, outpath):
    if len(np.shape(columns)) == 1:
        print("Only 1 col provided")
        print(np.shape(columns))
        columns = [columns]
    print("col shape: ", np.shape(columns))
    table = np.column_stack(columns)
    print(table.shape)
    np.savetxt(outpath, table, header=header, comments="")


def round_to_one_significant(x):
    """
    Rounds to one significant digit. 
    """
    return round(x, -int(np.floor(np.log10(abs(x)))))


def small_angle_corr_analytical(theta_max):
    """From semi guess aided analytical computation"""
    return (1 - np.cos(theta_max)) / 2


def small_angle_corr_litteratur(theta_max):
    """Approximation(?) from litterature, valid for small angles"""
    return (theta_max**2) / 6


def stepsize_isotropic(D, v):
    return 3 * D / v


def stepsize_analytical(theta_max, D, v):
    return stepsize_isotropic(D, v) * small_angle_corr_analytical(theta_max)


def power_law(x, a, b):
    return a * np.power(x, b)


def piecewise_poly(x, bp, coeffs):
    """
    Callable piecewise polynomial from breakpoints and coefficients

    parameters:
        x: real number, function argument
        bp: 1d array, breakpoints
        coeffs: ndarray, coefficients for polynomials on each interval
    """
    if x < bp[0] or x > bp[-1]:
        raise Exception(f"Error, argument x={x} out of range")
    for (i, b) in enumerate(bp[1:]):
        if x <= b:
            k = 3
            return sum(coeffs[m, i] * (x - bp[i]) ** (k - m) for m in range(k + 1))
    raise Exception(f"Unknown error, possibly invalid argument")


piecewise_poly_vectorized = np.vectorize(piecewise_poly)
piecewise_poly_vectorized.excluded.add(1)
piecewise_poly_vectorized.excluded.add(2)


def curve_fit_and_plot(verbose=False, include_analyt=False, plot=True, savedat=False):
    v = 0.99558849759333801  # Proton of energy 1GeV
    R_L = 3.5230000000000007e-005  # Larmor radius (isotropic stepsize)
    D = R_L / 3
    # print(f"D theory: {D}")

    theta_max_data_arr, D_coeff_arr = np.genfromtxt(in_file, delimiter=",").T
    # print(theta_max_data_arr)
    # print(D_coeff_arr)
    theta_max_data_arr *= np.pi
    theta_max_arr = np.linspace(theta_max_data_arr[0], theta_max_data_arr[-1], 100)

    corr_arr = (R_L * v) / (3 * D_coeff_arr)
    # D_target = D#1.1687289275261758e-05 # Replace with theoretical
    # print("Target D_coeff: ", D_target)

    # Fitting power law
    pars, cov = curve_fit(
        f=power_law,
        xdata=theta_max_data_arr[0:-10],
        ydata=corr_arr[0:-10],
        p0=[0, 0],
        bounds=(-np.inf, np.inf),
    )
    a, b = pars

    # Interpolation of stepsize (linear and cubic splines)
    l_ls = interp1d(theta_max_data_arr, corr_arr)  # linear interpolation
    l_cs = CubicSpline(theta_max_data_arr, corr_arr)  # Cubic splines

    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.set_xlabel(r"$\theta_{\mathrm{max}}$")
        ax1.set_ylabel("$\lambda / \lambda_0$")
        ax1.plot(theta_max_data_arr, corr_arr, "+", label="Data points", markersize=5)
        # ax1.plot(theta_max_arr, l_ls(theta_max_arr), label="Linear interpolation")
        ax1.plot(theta_max_arr, l_cs(theta_max_arr), label="Cubic Spline")
        ax1.plot(theta_max_arr, power_law(theta_max_arr, a, b), label="Power law fit")
        if include_analyt:
            ax1.plot(
                theta_max_arr,
                small_angle_corr_analytical(theta_max_arr),
                label="analytical attempt",
            )
        ax1.plot(
            theta_max_arr, small_angle_corr_litteratur(theta_max_arr), label="litterature"
        )
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.legend()
        ax2.set_xlabel(r"$\theta_{\mathrm{max}}$")
        # ax2.set_ylabel("$\lambda / \lambda_0$")
        ax2.plot(theta_max_data_arr, corr_arr, "+", label="Data points", markersize=5)
        # ax2.plot(theta_max_arr, l_ls(theta_max_arr), label="Linear interpolation")
        ax2.plot(theta_max_arr, l_cs(theta_max_arr), label="Cubic Spline")
        ax2.plot(theta_max_arr, power_law(theta_max_arr, a, b), label="Power law fit")
        if include_analyt:
            ax2.plot(
                theta_max_arr,
                small_angle_corr_analytical(theta_max_arr),
                label="analytical attempt",
            )
        ax2.plot(
            theta_max_arr, small_angle_corr_litteratur(theta_max_arr), label="litterature"
        )
        ax2.legend()
        if __name__ == "__main__":
            plt.show()
    if savedat:
        columns = [ 
            theta_max_arr, 
            l_cs(theta_max_arr),
            power_law(theta_max_arr, a, b),
            small_angle_corr_analytical(theta_max_arr),
            small_angle_corr_litteratur(theta_max_arr),
        ]
        header = "theta cubicspline powerlaw analyt litterature"
        save_table(columns, header, f"{OUT_DIR}small_angle_stepsize_correction.dat")

        columns_dp = [
            theta_max_data_arr, 
            corr_arr,
        ]
        header_dp = "theta corr"
        save_table(columns_dp, header_dp, f"{OUT_DIR}small_angle_stepsize_correction_DP.dat")

    # Print power law parameters
    # print("######################")
    # print("Power law parameters: ")
    # print("Power Law: f(x) = a*x^b")
    # print("a: ", a)
    # print("b: ", b)
    # print("######################")
    # print()

    ## Print l_cs breakpoints and coeffs
    # print("######################")
    # print("l_cs breakpoints: ", l_cs.x)
    # print("l_cs breakpoints shape: ", l_cs.x.shape)
    # print("l_cs coeffs: ", l_cs.c)
    # print("l_cs coeffs shape: ", l_cs.c.shape)
    # print("######################")


if __name__ == "__main__":
    curve_fit_and_plot(verbose=True, plot=True, savedat=True)
