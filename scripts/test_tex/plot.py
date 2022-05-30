import numpy as np
import matplotlib.pyplot as plt


def plot1():
    """Creates a simple plot
    :returns: TODO

    """
    theta = 22
    x = np.linspace(0, 2*np.pi, 100)
    y = np.sin(x)
    plt.plot(x, y, label=fr"$\theta = {1+1}$")
    plt.legend()
