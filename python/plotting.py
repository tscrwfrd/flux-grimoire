import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np


def plot_fct_dam_1d():
    """ """
    data = np.genfromtxt("./data/dam_break.csv", dtype=np.float64, delimiter=",")
    nrows, ncols = data.shape
    nrows = int(nrows / 2)
    xa = range(ncols)

    fig, [ax1, ax2] = plt.subplots(
        2, 1, sharex=True, layout="constrained", dpi=200, figsize=(8, 6)
    )

    (line1,) = ax1.plot([], [], "b-", linewidth=3)
    (line2,) = ax2.plot([], [], "r-", linewidth=3)

    ax1.set_xlim(0, ncols)
    ax1.set_ylim(0, 15)

    ax2.set_ylim(np.min(data[1::2]), np.max(data[1::2]))
    ax2.set_ylim(0, 30)

    ax1.set_ylabel("rho")
    ax2.set_ylabel("velocity")

    def update(frame):
        line1.set_data(xa, data[frame * 2])
        line2.set_data(xa, data[frame * 2 + 1])
        return line1, line2

    anim = FuncAnimation(fig, update, frames=int(len(data) / 2), interval=20)
    anim.save("./data/fct_dam.gif")
    plt.close()

def plot_square_wave():
    """ """
    data = np.genfromtxt("./data/square_wave.csv", dtype=np.float64, delimiter=",")
    nrows, ncols = data.shape
    nrows = int(nrows / 2)
    xa = range(ncols)

    fig, [ax1, ax2] = plt.subplots(
        2, 1, sharex=True, layout="constrained", dpi=200, figsize=(8, 6)
    )

    (line1,) = ax1.plot([], [], "b-", linewidth=3)
    (line2,) = ax2.plot([], [], "r-", linewidth=3)

    ax1.set_xlim(0, ncols)
    ax1.set_ylim(0, 300.0)

    ax2.set_ylim(np.min(data[1::2]), np.max(data[1::2]))
    ax2.set_ylim(0, 1100.0)

    ax1.set_ylabel("rho")
    ax2.set_ylabel("velocity")

    def update(frame):
        line1.set_data(xa, data[frame * 2])
        line2.set_data(xa, data[frame * 2 + 1])
        return line1, line2

    anim = FuncAnimation(fig, update, frames=int(len(data) / 2), interval=20)
    anim.save("./data/square_wave.gif")
    plt.close()

if __name__ == "__main__":
    plot_square_wave()
    plot_fct_dam_1d()
