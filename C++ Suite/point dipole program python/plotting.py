"""
plotting.py
===========

Visualization utilities for Dyson orbitals and wavefunctions, translated
from MATLAB figplot.m.

Functions
---------
figplot(X, Y, Z, F, axval, iso, no, LR, ptype)
    Plot isosurface(s) of a 3-D function using matplotlib.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

__all__ = ["figplot"]


def figplot(
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    F: np.ndarray,
    axval: float,
    iso: float,
    no: int = 1,
    LR: int = 1,
    ptype: int = 1,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Plot positive and negative isosurfaces of a 3-D function.

    Parameters
    ----------
    X, Y, Z : ndarray (npts, npts, npts)
        Cartesian coordinate meshes.
    F : ndarray (npts, npts, npts)
        Function values (e.g., Dyson orbital).
    axval : float
        Axis limit (+/- axval for x, y, z).
    iso : float
        Isovalue for the isosurface.
    no : int, optional
        Number of subplots in the figure (default 1).
    LR : int, optional
        Subplot index (1-indexed as in MATLAB, default 1).
    ptype : int, optional
        Plot type flag; 1 = draw isosurface (default).
    ax : Axes, optional
        Existing 3-D axes to draw on.  If None, a new figure is created.

    Returns
    -------
    ax : Axes
        The matplotlib 3-D axes with the plot.

    Notes
    -----
    MATLAB source (figplot.m) uses MATLAB's ``isosurface`` and ``patch``.
    Here we use ``skimage.measure.marching_cubes`` to generate mesh
    vertices/faces and ``Poly3DCollection`` to render them.
    """
    if ptype != 1:
        raise NotImplementedError("Only ptype=1 (isosurface) is supported")

    # Create axes if not provided
    if ax is None:
        fig = plt.figure(figsize=(6 * no, 6))
        ax = fig.add_subplot(1, no, LR, projection="3d")

    # Extract grid spacing
    dx = X[1, 0, 0] - X[0, 0, 0]
    dy = Y[0, 1, 0] - Y[0, 0, 0]
    dz = Z[0, 0, 1] - Z[0, 0, 0]

    # Helper to add an isosurface
    def _add_iso(level: float, color: str, alpha: float = 0.3):
        try:
            verts, faces, _, _ = measure.marching_cubes(
                F, level=level, spacing=(dx, dy, dz)
            )
        except ValueError:
            return  # no surface at this level
        # Shift verts to match X, Y, Z origin
        verts[:, 0] += X.min()
        verts[:, 1] += Y.min()
        verts[:, 2] += Z.min()
        mesh = Poly3DCollection(verts[faces], alpha=alpha, edgecolor=color, facecolor="none")
        ax.add_collection3d(mesh)

    _add_iso(+iso, "red")
    _add_iso(-iso, "blue")

    # Axis settings
    ax.set_xlim(-axval, axval)
    ax.set_ylim(-axval, axval)
    ax.set_zlim(-axval, axval)
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()

    # Draw coordinate axes through origin
    ax.plot([0, 0], [0, 0], [-8, 8], "k.-")
    ax.plot([-4, 4], [0, 0], [0, 0], "k.-")
    ax.plot([0, 0], [-4, 4], [0, 0], "k.-")

    ax.view_init(elev=15, azim=30)
    return ax
