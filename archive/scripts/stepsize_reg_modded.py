import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit


def small_angle_corr_analytical(theta_max):
    """ From semi guess aided analytical computation """
    return (1-np.cos(theta_max))/2


def stepsize_isotropic(D, v):
    return 3*D/v


def stepsize_analytical(theta_max, D, v):
    return stepsize_isotropic(D, v) * small_angle_corr_analytical(theta_max)


def power_law(x, a, b):
    return a*np.power(x, b)


def piecewise_poly(x, bp, coeffs):
    """ 
    Callable piecewise polynomial from breakpoints and coefficients 

    parameters:
        x: real number, function argument
        bp: 1d array, breakpoints
        coeffs: ndarray, coefficients for polynomials on each interval
    """

#    print("LARGEST BP")
#    print(max(bp))
#
#    print("piecewise_poly paramaters: ")
#    print("####################")
#    print("x: ", x)
#    print()
#    print("bp (breakpoints): \n")
#    print(bp)
#    print("")
#    print("coeffs (coefficients): \n")
#    print(coeffs)
#    print("####################")

    if x < bp[0] or x > bp[-1]:
        raise Exception(f"Error, argument x={x} out of range")
    for (i, b) in enumerate(bp[1:]):
        if x <= b:
            k = 3
            if x > bp[-2]:
                print("x>0.9pi BBY:")
                print("coeffs")
                print(coeffs[:, i])
                print("bp")
                print(bp[i])
            return sum(coeffs[m, i] * (x - bp[i])**(k-m) for m in range(k+1))
    raise Exception(f"Unknown error, possibly invalid argument")


piecewise_poly_vectorized = np.vectorize(piecewise_poly)
piecewise_poly_vectorized.excluded.add(1)
piecewise_poly_vectorized.excluded.add(2)


if __name__ == "__main__":
    v = 0.99558849759333801 # Proton of energy 1GeV
    R_L = 3.5230000000000007E-005 # Larmor radius (isotropic stepsize)
    D = R_L/3 # Target diffusion coeficient

    # Data
    theta_max_data_arr, D_coeff_arr = np.genfromtxt("rw_dcoeff.csv", delimiter=",").T
    theta_max_data_arr *= np.pi
    # 'theory' values
    theta_max_arr = np.linspace(0.001, theta_max_data_arr[-1], 100)
    theta_max_arr *= np.pi

    corr_arr = R_L*v/(3*D_coeff_arr)
    #D_target = 1.1750165596133696e-05
    #print("Target D_coeff: ", D_target)

    # Fitting power law (seems valid for small enough angles):
    pars, cov = curve_fit(f=power_law, xdata=theta_max_data_arr, ydata=corr_arr, p0=[0, 0], bounds=(-np.inf, np.inf))
    stdevs = np.sqrt(np.diag(cov))
    a, b = pars
    print("##################")
    print("Powerlaw: f(x) = a*x^b")
    print("a, b: ", a, b)
    print("##################")

    # Simple slope calculation to compare
    ln_theta_max_data_arr = np.log(theta_max_data_arr)
    ln_corr_arr = np.log(corr_arr)
    power = (ln_corr_arr[-1]-ln_corr_arr[0])/(ln_theta_max_data_arr[-1]-ln_theta_max_data_arr[0])
    print(f"power_cf - power_slope: {b - power}")

    plt.title("Small angle stepsize correction factor -- powerlaw fit (log scale)")
    plt.plot(theta_max_data_arr, corr_arr, "o", label="Data points")
    plt.plot(theta_max_data_arr, power_law(theta_max_data_arr, a, b), 'x', label="Fitted disc")
    plt.plot(theta_max_arr, power_law(theta_max_arr, a, b), label="Fitted cont.")
    plt.plot(theta_max_arr, small_angle_corr_analytical(theta_max_arr), label="analytical attempt")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.xlabel("theta max (largest scattering angle)")
    plt.ylabel("stepsize correction f(theta_max)")
    plt.show()

    plt.title("Small angle stepsize correction factor -- powerlaw fit (normal scale)")
    plt.plot(theta_max_data_arr, corr_arr, "o", label="Data points")
    plt.plot(theta_max_data_arr, power_law(theta_max_data_arr, a, b), 'x', label="Fitted disc")
    plt.plot(theta_max_arr, power_law(theta_max_arr, a, b), label="Fitted cont.")
    plt.plot(theta_max_arr, small_angle_corr_analytical(theta_max_arr), label="analytical attempt")
    plt.legend()
    plt.xlabel("theta max (largest scattering angle)")
    plt.ylabel("stepsize correction f(theta_max)")
    plt.show()
