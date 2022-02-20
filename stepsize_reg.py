import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline


def stepsize_analytical(theta_max, D, v):
    return (3*D/v) * (1 - np.cos(theta_max))/2


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
    D = R_L/3
    theta_max_data_arr, D_prime_coeff_arr = np.genfromtxt("rw_dcoeff.csv", delimiter=",").T
    theta_max_arr = np.linspace(theta_max_data_arr[0], theta_max_data_arr[-1], 100)
    theta_max_data_arr *= np.pi
    theta_max_arr *= np.pi
    K_arr = 3*D_prime_coeff_arr/(R_L*v)
    D_prime_target = D_prime_coeff_arr[-1]
    print("Target D_coeff: ", D_prime_target)

    #plt.plot(theta_max_arr, D_prime_coeff_arr)
    #plt.show()

#    # Exponential fit K(theta_max_arr)#
#    ln_K_arr = np.log(K_arr)
#    theta_max_data_inv_arr = 1/theta_max_data_arr
#    K_ln = np.polyfit(theta_max_data_inv_arr, ln_K_arr, 1)
#    A = np.exp(K_ln[1])
#    B = K_ln[0]
#
#    # spline interpolation of K
#    K_ls = interp1d(theta_max_data_arr, K_arr)
#    K_cs = CubicSpline(theta_max_data_arr, K_arr)
#    #cs = interp1d(theta_max_data_arr, K_arr, kind="cubic")
#
#    # Stepsize estimate via K
#    l_K_ls = R_L/K_ls(theta_max_arr)
#    l_K_cs = R_L/K_cs(theta_max_arr)

    # Plot
    #plt.title("K(theta_max)")
    #plt.plot(theta_max_data_arr, K_arr, ".", label="Data")
    #plt.plot(theta_max_arr, A*np.exp(B/theta_max_arr), label="exp fit")
    #plt.plot(theta_max_arr, K_ls(theta_max_arr), label="LS")
    #plt.plot(theta_max_arr, K_cs(theta_max_arr), label="CS")
    #plt.legend()
    #plt.show()

    ######################################################################

    # Direct fit/interpolation of stepsize
    l_ls = interp1d(theta_max_data_arr, R_L/K_arr) # linear interpolation
    l_cs = CubicSpline(theta_max_data_arr, R_L/K_arr) # Cubic splines
    p2_coeffs = np.polyfit(theta_max_data_arr, R_L/K_arr, 2) # polyfit (deg 2)
    p3_coeffs = np.polyfit(theta_max_data_arr, R_L/K_arr, 3) # polyfit (deg 3)
    l_poly2 = lambda x: p2_coeffs[0]*x**2 + p2_coeffs[1]*x + p2_coeffs[2]
    l_poly3 = lambda x: p3_coeffs[0]*x**3 + p3_coeffs[1]*x**2 + p3_coeffs[2]*x + p3_coeffs[3]

    plt.title("stepsize asf. theta_max")
    plt.xlabel("theta_max")
    plt.ylabel("stepsize")
    plt.plot(theta_max_data_arr, R_L/K_arr, "o", label="DP")
    #plt.plot(theta_max_arr, l_ls(theta_max_arr), label="Linear interpolation")
    plt.plot(theta_max_arr, l_cs(theta_max_arr), label="Cubic Spline")
    #plt.plot(theta_max_arr, stepsize_analytical(theta_max_arr, D, v), label='analytical attempt')
    #plt.plot(theta_max_arr, l_poly2(theta_max_arr), label="2deg polyfit")
    #plt.plot(theta_max_arr, l_poly3(theta_max_arr), label="3deg polyfit")
    plt.plot(theta_max_arr, stepsize_analytical(theta_max_arr, D, v)/stepsize_analytical(theta_max_arr, D, v)[-1], label='analytical attempt')
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.show()

#    plt.plot(theta_max_arr, piecewise_poly_vectorized(theta_max_arr, l_cs.x, l_cs.c), label="qubic spline")
##    plt.plot(theta_max_arr, stepsize_analytical(theta_max_arr, D, v)/stepsize_analytical(theta_max_arr, D, v)[-1], label='analytical attempt')
#    plt.legend()
#    plt.show()


    plt.title("Small angle reduction factor")
    plt.plot(theta_max_arr, stepsize_analytical(theta_max_arr, D, v)/stepsize_analytical(theta_max_arr, D, v)[-1], label='analytical attempt')
    plt.plot(theta_max_arr, piecewise_poly_vectorized(theta_max_arr, l_cs.x, l_cs.c)/piecewise_poly_vectorized(theta_max_arr, l_cs.x, l_cs.c)[-1], label="qubic spline")
    plt.legend()
    plt.show()

    #print("FACT: ", piecewise_poly_vectorized(theta_max_arr, l_cs.x, l_cs.c)[-1])

#    # Print l_cs breakpoints and coeffs
#    print("l_cs breakpoints: ", l_cs.x)
#    print("l_cs breakpoints shape: ", l_cs.x.shape)
#    print("l_cs coeffs: ", l_cs.c)
#    print("l_cs coeffs shape: ", l_cs.c.shape)
#    # Write function/subroutine in Fortran to construct callable from list of breakpoints and coeffs


    #print(piecewise_poly_vectorized(0.005*np.pi, l_cs.x, l_cs.c))
    #print(piecewise_poly_vectorized(0.1*np.pi, l_cs.x, l_cs.c))
    #print(piecewise_poly_vectorized(1.0*np.pi, l_cs.x, l_cs.c))
