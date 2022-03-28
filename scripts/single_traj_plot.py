import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    frac = 0.1
    data = np.genfromtxt(f"./data/out-shock{frac:.1f}")
    x = data.T[0]
    y = data.T[1]
    z = data.T[2]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    plt.title(f"theta-max = {frac} * pi")
    ax.plot(x, y, z)
    ax.plot(x, y, z, ".")
    plt.show()
