#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as mpl
import numpy as np
import pymap3d as pm

p = argparse.ArgumentParser()
p.add_argument("alt_m", help="altitude [meters]", type=float, default=0.0, nargs="?")
args = p.parse_args()

lat, lon = np.meshgrid(np.arange(-90, 90, 0.1), np.arange(-180, 180, 0.2))

x, y, z = pm.geodetic2ecef(lat, lon, args.alt_m)


def panel(ax, val, name: str, cmap: str = None):
    hi = ax.pcolormesh(lon, lat, val, cmap=cmap)
    ax.set_title(name)
    fg.colorbar(hi, ax=ax).set_label(name + " [m]")
    ax.set_xlabel("longitude [deg]")


fg = mpl.figure(figsize=(16, 5))
axs = fg.subplots(1, 3, sharey=True)
fg.suptitle("geodetic2ecef")

panel(axs[0], x, "x", "bwr")
panel(axs[1], y, "y", "bwr")
panel(axs[2], z, "z", "bwr")

axs[0].set_ylabel("latitude [deg]")


fg = mpl.figure(figsize=(16, 5))
axs = fg.subplots(1, 3, sharey=True)
fg.suptitle(r"|$\nabla$ geodetic2ecef|")


panel(axs[0], np.hypot(*np.gradient(x)), r"|$\nabla$ x|")
panel(axs[1], np.hypot(*np.gradient(y)), r"|$\nabla$ y|")
panel(axs[2], np.hypot(*np.gradient(z)), r"|$\nabla$ z|")

axs[0].set_ylabel("latitude [deg]")

mpl.show()
