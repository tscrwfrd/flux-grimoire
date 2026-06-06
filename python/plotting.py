import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np


def plot_sod_shock_1d_lw():
    """ """
    data = np.genfromtxt("./data/sod_shock_lw.csv", dtype=np.float64, delimiter=",")
    nrows, ncols = data.shape
    nrows = int(nrows / 3)
    xa = range(ncols)

    fig, [ax1, ax2, ax3, ax4] = plt.subplots(
        4, 1, sharex=True, layout="constrained", dpi=200, figsize=(8, 8)
    )

    (line1,) = ax1.plot([], [], "b-", linewidth=3)
    (line2,) = ax2.plot([], [], "r-", linewidth=3)
    (line3,) = ax3.plot([], [], "g-", linewidth=3)
    (line4,) = ax4.plot([], [], "m-", linewidth=3)

    ax1.set_xlim(0, ncols)
    ax1.set_ylim(0, 1.2)

    ax2.set_ylim(np.min(data[1::2]), np.max(data[1::2]))
    ax2.set_ylim(0, 2)

    ax3.set_ylim(0, 1.5)
    ax4.set_ylim(0, 3.0)

    ax1.set_ylabel("rho")
    ax2.set_ylabel("velocity")
    ax3.set_ylabel("pressure")
    ax4.set_ylabel("energy")

    def update(frame):
        rho = data[frame * 3]
        vel = data[frame * 3 + 1]
        prs = data[frame * 3 + 2]
        energy = prs / (1.4 - 1.0) + 0.5 * rho * vel**2
        line1.set_data(xa, rho)
        line2.set_data(xa, vel)
        line3.set_data(xa, prs)
        line4.set_data(xa, energy)
        return line1, line2, line3, line4

    anim = FuncAnimation(fig, update, frames=int(len(data) / 3), interval=100)
    anim.save("./data/sod_shock_lw.gif")
    plt.close()

def plot_sod_shock_1d_roe():
    """ """
    data = np.genfromtxt("./data/sod_shock_roe.csv", dtype=np.float64, delimiter=",")
    nrows, ncols = data.shape
    nrows = int(nrows / 3)
    xa = range(ncols)

    fig, [ax1, ax2, ax3, ax4] = plt.subplots(
        4, 1, sharex=True, layout="constrained", dpi=200, figsize=(8, 8)
    )

    (line1,) = ax1.plot([], [], "b-", linewidth=3)
    (line2,) = ax2.plot([], [], "r-", linewidth=3)
    (line3,) = ax3.plot([], [], "g-", linewidth=3)
    (line4,) = ax4.plot([], [], "m-", linewidth=3)

    ax1.set_xlim(0, ncols)
    ax1.set_ylim(0, 1.2)

    ax2.set_ylim(np.min(data[1::2]), np.max(data[1::2]))
    ax2.set_ylim(0, 2)

    ax3.set_ylim(0, 1.5)
    ax4.set_ylim(0, 3.0)

    ax1.set_ylabel("rho")
    ax2.set_ylabel("velocity")
    ax3.set_ylabel("pressure")
    ax4.set_ylabel("energy")

    def update(frame):
        rho = data[frame * 3]
        vel = data[frame * 3 + 1]
        prs = data[frame * 3 + 2]
        energy = prs / (1.4 - 1.0) + 0.5 * rho * vel**2
        line1.set_data(xa, rho)
        line2.set_data(xa, vel)
        line3.set_data(xa, prs)
        line4.set_data(xa, energy)
        return line1, line2, line3, line4

    anim = FuncAnimation(fig, update, frames=int(len(data) / 3), interval=100)
    anim.save("./data/sod_shock_roe.gif")
    plt.close()

def plot_sod_shock_1d_weno3():
    """ """
    data = np.genfromtxt("./data/sod_shock_weno3.csv", dtype=np.float64, delimiter=",")
    nrows, ncols = data.shape
    nrows = int(nrows / 3)
    xa = range(ncols)

    fig, [ax1, ax2, ax3, ax4] = plt.subplots(
        4, 1, sharex=True, layout="constrained", dpi=200, figsize=(8, 8)
    )

    (line1,) = ax1.plot([], [], "b-", linewidth=3)
    (line2,) = ax2.plot([], [], "r-", linewidth=3)
    (line3,) = ax3.plot([], [], "g-", linewidth=3)
    (line4,) = ax4.plot([], [], "m-", linewidth=3)

    ax1.set_xlim(0, ncols)
    ax1.set_ylim(0, 1.2)

    ax2.set_ylim(np.min(data[1::2]), np.max(data[1::2]))
    ax2.set_ylim(0, 2)

    ax3.set_ylim(0, 1.5)
    ax4.set_ylim(0, 3.0)

    ax1.set_ylabel("rho")
    ax2.set_ylabel("velocity")
    ax3.set_ylabel("pressure")
    ax4.set_ylabel("energy")

    def update(frame):
        rho = data[frame * 3]
        vel = data[frame * 3 + 1]
        prs = data[frame * 3 + 2]
        energy = prs / (1.4 - 1.0) + 0.5 * rho * vel**2
        line1.set_data(xa, rho)
        line2.set_data(xa, vel)
        line3.set_data(xa, prs)
        line4.set_data(xa, energy)
        return line1, line2, line3, line4

    anim = FuncAnimation(fig, update, frames=int(len(data) / 3), interval=100)
    anim.save("./data/sod_shock_weno3.gif")
    plt.close()

def plot_sod_shock_1d_weno5():
    """ """
    data = np.genfromtxt("./data/sod_shock_weno5.csv", dtype=np.float64, delimiter=",")
    nrows, ncols = data.shape
    nrows = int(nrows / 3)
    xa = range(ncols)

    fig, [ax1, ax2, ax3, ax4] = plt.subplots(
        4, 1, sharex=True, layout="constrained", dpi=200, figsize=(8, 8)
    )

    (line1,) = ax1.plot([], [], "b-", linewidth=3)
    (line2,) = ax2.plot([], [], "r-", linewidth=3)
    (line3,) = ax3.plot([], [], "g-", linewidth=3)
    (line4,) = ax4.plot([], [], "m-", linewidth=3)

    ax1.set_xlim(0, ncols)
    ax1.set_ylim(0, 1.2)

    ax2.set_ylim(np.min(data[1::2]), np.max(data[1::2]))
    ax2.set_ylim(0, 2)

    ax3.set_ylim(0, 1.5)
    ax4.set_ylim(0, 3.0)

    ax1.set_ylabel("rho")
    ax2.set_ylabel("velocity")
    ax3.set_ylabel("pressure")
    ax4.set_ylabel("energy")

    def update(frame):
        rho = data[frame * 3]
        vel = data[frame * 3 + 1]
        prs = data[frame * 3 + 2]
        energy = prs / (1.4 - 1.0) + 0.5 * rho * vel**2
        line1.set_data(xa, rho)
        line2.set_data(xa, vel)
        line3.set_data(xa, prs)
        line4.set_data(xa, energy)
        return line1, line2, line3, line4

    anim = FuncAnimation(fig, update, frames=int(len(data) / 3), interval=100)
    anim.save("./data/sod_shock_weno5.gif")
    plt.close()

def plot_fct_dam_1d():
    """ """
    data = np.genfromtxt("./data/dam_break.csv", dtype=np.float64, delimiter=",")
    nrows, ncols = data.shape
    nrows = int(nrows / 2)
    xa = range(ncols)

    fig, [ax1, ax2, ax3] = plt.subplots(
        3, 1, sharex=True, layout="constrained", dpi=200, figsize=(8, 7)
    )

    (line1,) = ax1.plot([], [], "b-", linewidth=3)
    (line2,) = ax2.plot([], [], "r-", linewidth=3)
    (line3,) = ax3.plot([], [], "m-", linewidth=3)

    ax1.set_xlim(0, ncols)
    ax1.set_ylim(0, 4)

    ax2.set_ylim(np.min(data[1::2]), np.max(data[1::2]))
    ax2.set_ylim(0, 5)

    ax3.set_ylim(0, 50)

    ax1.set_ylabel("rho")
    ax2.set_ylabel("velocity")
    ax3.set_ylabel("energy")

    g = 9.8

    def update(frame):
        hgt = data[frame * 2]
        vel = data[frame * 2 + 1]
        energy = 0.5 * hgt * vel**2 + 0.5 * g * hgt**2
        line1.set_data(xa, hgt)
        line2.set_data(xa, vel)
        line3.set_data(xa, energy)
        return line1, line2, line3

    anim = FuncAnimation(fig, update, frames=int(len(data) / 2), interval=100)
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

def _animate_2d_field(csv_path, gif_path, vmin=0.0, vmax=1.0):
    """Read a 2D-per-snapshot CSV and render it as an animated heatmap."""
    data = np.genfromtxt(csv_path, dtype=np.float64, delimiter=",")
    total_rows, nx = data.shape
    ny = nx
    nframes = total_rows // ny
    frames = data.reshape(nframes, ny, nx)

    fig, ax = plt.subplots(layout="constrained", dpi=150, figsize=(6, 6))
    im = ax.imshow(
        frames[0],
        origin="lower",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    fig.colorbar(im, ax=ax, shrink=0.8, label="rho")

    def update(frame):
        im.set_data(frames[frame])
        ax.set_title(f"frame {frame}/{nframes - 1}")
        return (im,)

    anim = FuncAnimation(fig, update, frames=nframes, interval=120)
    anim.save(gif_path)
    plt.close()


def plot_slotted_cylinder():
    """Animate Zalesak's 2D slotted-cylinder rotation test (FCT result).

    CSV layout: ny consecutive rows of nx values per snapshot, square grid
    so ny is recovered as nx (rows-per-frame).
    """
    _animate_2d_field(
        "./data/slotted_cylinder.csv", "./data/slotted_cylinder.gif"
    )

if __name__ == "__main__":
    plot_square_wave()
    plot_fct_dam_1d()
    plot_sod_shock_1d_lw()
    plot_sod_shock_1d_roe()
    plot_sod_shock_1d_weno3()
    plot_sod_shock_1d_weno5()
    plot_slotted_cylinder()
