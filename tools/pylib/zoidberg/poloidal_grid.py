"""Routines for generating structured meshes on poloidal domains

Classes
-------
`RectangularPoloidalGrid`
    Simple rectangles in R-Z
`StructuredPoloidalGrid`
    Curvilinear structured grids in R-Z

Functions
---------
`grid_annulus`
    Create a `StructuredPoloidalGrid` from inner and outer RZLine
    objects using a simple algorithm
`grid_elliptic`
    Create a `StructuredPoloidalGrid` from inner and outer RZLine
    objects using elliptic meshing method

"""

import numpy as np
from numpy import pi, linspace, zeros
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import cKDTree as KDTree

import warnings

try:
    import matplotlib.pyplot as plt

    plotting_available = True
except ImportError:
    warnings.warn("Couldn't import matplotlib, plotting not available.")
    plotting_available = False

try:
    from . import rzline
except:
    # Python 2
    import rzline


class PoloidalGrid(object):
    """Represents a poloidal grid

    Note: Here the 2D plane (R,Z) is labelled by (x,z) indices

    Attributes
    ----------
    nx, nz : int
        Number of points in x and z
    R : ndarray
        2D Numpy array of R coordinates
    Z : ndarray
        2D Numpy array of Z coordinates

    """

    def plot(self, axis=None, show=True):
        """Plot grid using matplotlib

        Parameters
        ----------
        axis : matplotlib axis, optional
            A matplotlib axis to plot on. By default a new figure
            is created
        show : bool, optional
            Calls plt.show() at the end

        Returns
        -------
        axis
            The matplotlib axis that was used

        """

        if not plotting_available:
            warnings.warn("matplotlib not available, unable to plot")
            return None

        if axis is None:
            fig = plt.figure()
            axis = fig.add_subplot(1, 1, 1)

        axis.plot(self.R, self.Z, "k-")
        axis.plot(self.R.T, self.Z.T, "k-")

        if show:
            plt.show()

        return axis


class RectangularPoloidalGrid(PoloidalGrid):
    """Represents a poloidal grid consisting of a rectangular domain

    Note: Here the 2D plane (R,Z) is labelled by (x,z) indices

    Attributes
    ----------
    nx, nz : int
        Number of points in x and z
    R : ndarray
        2D Numpy array of R coordinates
    Z : ndarray
        2D Numpy array of Z coordinates

    Parameters
    ----------
    nx : int
        Number of points in major radius (including boundaries)
    nz : int
        Number of points in height (including boundaries)
    Lx : float
        Radial domain size  [m]
    Lz : float
        Vertical domain size [m]
    Rcentre : float, optional
        Coordinate at the middle of the domain
    Zcentre : float, optional
        Coordinate at the middle of the domain
    MXG : int, optional
        Number of guard cells in X. The boundary is put half-way
        between the guard cell and the domain

    """

    def __init__(self, nx, nz, Lx, Lz, Rcentre=0.0, Zcentre=0.0, MXG=2):

        self.nx = nx
        self.nz = nz

        self.Lx = Lx
        self.Lz = Lz

        self.Rcentre = Rcentre
        self.Zcentre = Zcentre

        # Some useful derived quantities
        # Note: index at the middle of the domain is (nx - 1)/2
        #       e.g. nx=5 :  0 1 | 2 | 3 4
        self.dR = self.Lx / (self.nx - 2 * MXG)
        self.dZ = self.Lz / self.nz
        self.Rmin = self.Rcentre - self.dR * (self.nx - 1.0) / 2.0
        self.Zmin = self.Zcentre - self.dZ * (self.nz - 1.0) / 2.0

        # Generate 2D arrays
        # Using getCoordinate to ensure consistency
        xind, zind = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
        self.R, self.Z = self.getCoordinate(xind, zind)

    def __repr__(self):
        return (
            "RectangularPoloidalGrid({0},{1},{2},{3},Rcentre={4},Zcentre={5})".format(
                self.nx, self.nz, self.Lx, self.Lz, self.Rcentre, self.Zcentre
            )
        )

    def getCoordinate(self, xind, zind, dx=0, dz=0):
        """Get coordinates (R,Z) at given (xind,zind) index

        Parameters
        ----------
        xind, zind : array_like
            Indices in X and Z. These should be the same shape
        dx : int, optional
            Order of x derivative
        dz : int, optional
            Order of z derivative

        Returns
        -------
        R, Z : (ndarray, ndarray)
            Locations of point or derivatives of R,Z with respect to
            indices if dx,dz != 0

        """

        # Convert to NumPy arrays if not already
        xind = np.asfarray(xind)
        zind = np.asfarray(zind)
        # Make sure dx and dz are integers
        dx = int(dx)
        dz = int(dz)

        assert xind.shape == zind.shape
        assert dx >= 0
        assert dz >= 0

        shape = xind.shape

        if dx + dz > 2:
            # Second derivatives and mixed derivatives all zero
            return np.zeros(shape), np.zeros(shape)

        if dx == 1:
            # dR/dx, dZ/dx
            return np.full(shape, self.dR), np.zeros(shape)
        elif dz == 1:
            # dR/dz, dZ/dz
            return np.zeros(shape), np.full(shape, self.dZ)
        # Return (R,Z) location
        return self.Rmin + xind * self.dR, self.Zmin + zind * self.dZ

    def findIndex(self, R, Z):
        """Finds the (x,z) index corresponding to the given (R,Z) coordinate

        Parameters
        ----------
        R, Z : array_like
            Locations to find indices for

        Returns
        -------
        x, z : (ndarray, ndarray)
            Index as a float, same shape as R,Z

        """

        # Make sure inputs are NumPy arrays
        R = np.asfarray(R)
        Z = np.asfarray(Z)

        # Check that they have the same shape
        assert R.shape == Z.shape

        xind = (R - self.Rmin) / self.dR
        zind = (Z - self.Zmin) / self.dZ

        # Note: These indices may be outside the domain,
        # but this is handled in BOUT++, and useful for periodic
        # domains.

        return xind, zind

    def metric(self):
        """Return the metric tensor, dx and dz

        For this rectangular grid the metric is the identity

        Returns
        -------
        dict
            Dictionary containing:
            - **dx, dz**: Grid spacing
            - **gxx, gxz, gzz**: Covariant components
            - **g_xx, g_xz, g_zz**: Contravariant components
        """
        return {
            "dx": self.dR,
            "dz": self.dZ,  # Grid spacing
            "gxx": 1.0,
            "g_xx": 1.0,
            "gxz": 0.0,
            "g_xz": 0.0,
            "gzz": 1.0,
            "g_zz": 1.0,
        }


class StructuredPoloidalGrid(PoloidalGrid):
    """Represents a structured poloidal grid in R-Z

    Attributes
    ----------
    nx, nz : int
        Number of points in x and z
    R : ndarray
        2D Numpy array of R coordinates
    Z : ndarray
        2D Numpy array of Z coordinates

    Parameters
    ----------
    R, Z : ndarray
        2D Numpy arrays of R,Z points

        .. note:: R,Z are not copied, so these arrays should not be
                  modified afterwards

    """

    def __init__(self, R, Z):

        assert R.shape == Z.shape

        self.R = R
        self.Z = Z

        # Create a KDTree for quick lookup of nearest points
        n = R.size
        data = np.concatenate((R.reshape((n, 1)), Z.reshape((n, 1))), axis=1)
        self.tree = KDTree(data)

        # Create splines for quick interpolation of coordinates
        nx, nz = R.shape

        self.nx = nx
        self.nz = nz

        xinds = np.arange(nx)
        zinds = np.arange(nz + 1)
        # Repeat the final point in y since periodic in y
        R_ext = np.concatenate((R, np.reshape(R[:, 0], (nx, 1))), axis=1)
        Z_ext = np.concatenate((Z, np.reshape(Z[:, 0], (nx, 1))), axis=1)

        self._spl_r = RectBivariateSpline(xinds, zinds, R_ext)
        self._spl_z = RectBivariateSpline(xinds, zinds, Z_ext)

    def __repr__(self):
        return "StructuredPoloidalGrid()"

    def getCoordinate(self, xind, zind, dx=0, dz=0):
        """Get coordinates (R, Z) at given (xind, zind) index

        Parameters
        ----------
        xind, zind : array_like
            Indices in X and Z. These should be the same shape
        dx : int, optional
            Order of x derivative
        dz : int, optional
            Order of z derivative

        Returns
        -------
        R, Z : (ndarray, ndarray)
            Locations of point or derivatives of R,Z with respect to
            indices if dx,dz != 0

        """
        nx, nz = self.R.shape
        if (np.amin(xind) < 0) or (np.amax(xind) > nx - 1):
            raise ValueError("x index out of range")

        # Periodic in y
        zind = np.remainder(zind, nz)

        R = self._spl_r(xind, zind, dx=dx, dy=dz, grid=False)
        Z = self._spl_z(xind, zind, dx=dx, dy=dz, grid=False)

        return R, Z

    def findIndex(self, R, Z, tol=1e-10, show=False):
        """Finds the (x, z) index corresponding to the given (R, Z) coordinate

        Parameters
        ----------
        R, Z : array_like
            Locations. Can be scalar or array, must be the same shape
        tol : float, optional
            Maximum tolerance on the square distance

        Returns
        -------
        x, z : (ndarray, ndarray)
            Index as a float, same shape as R, Z

        """

        # Make sure inputs are NumPy arrays
        R = np.asfarray(R)
        Z = np.asfarray(Z)

        # Check that they have the same shape
        assert R.shape == Z.shape

        input_shape = R.shape  # So output has same shape as input

        # Get distance and index into flattened data
        # Note ind can be an integer, or an array of ints
        # with the same number of elements as the input (R,Z) arrays
        n = R.size
        position = np.concatenate((R.reshape((n, 1)), Z.reshape((n, 1))), axis=1)

        R = R.reshape((n,))
        Z = Z.reshape((n,))

        dists, ind = self.tree.query(position)

        # Calculate (x,y) index
        nx, nz = self.R.shape
        xind = np.floor_divide(ind, nz)
        zind = ind - xind * nz

        # Convert indices to float
        xind = np.asfarray(xind)
        zind = np.asfarray(zind)

        # Create a mask for the positions
        mask = np.ones(xind.shape)
        mask[
            np.logical_or((xind < 0.5), (xind > (nx - 1.5)))
        ] = 0.0  # Set to zero if near the boundary

        if show and plotting_available:
            plt.plot(self.R, self.Z, ".")
            plt.plot(R, Z, "x")

        while True:
            # Use Newton iteration to find the index
            # dR, dZ are the distance away from the desired point
            Rpos, Zpos = self.getCoordinate(xind, zind)
            if show and plotting_available:
                plt.plot(Rpos, Zpos, "o")
            dR = Rpos - R
            dZ = Zpos - Z

            # Check if close enough
            # Note: only check the points which are not in the boundary
            if np.amax(mask * (dR ** 2 + dZ ** 2)) < tol:
                break

            # Calculate derivatives
            dRdx, dZdx = self.getCoordinate(xind, zind, dx=1)
            dRdz, dZdz = self.getCoordinate(xind, zind, dz=1)

            # Invert 2x2 matrix to get change in coordinates
            #
            # (x) -=  ( dR/dx   dR/dz )^-1  (dR)
            # (y)     ( dZ/dx   dZ/dz )     (dz)
            #
            #
            # (x) -=  ( dZ/dz  -dR/dz ) (dR)
            # (y)     (-dZ/dx   dR/dx ) (dZ) / (dR/dx*dZ/dy - dR/dy*dZ/dx)
            determinant = dRdx * dZdz - dRdz * dZdx

            xind -= mask * ((dZdz * dR - dRdz * dZ) / determinant)
            zind -= mask * ((dRdx * dZ - dZdx * dR) / determinant)

            # Re-check for boundary
            in_boundary = xind < 0.5
            mask[in_boundary] = 0.0  # Set to zero if near the boundary
            xind[in_boundary] = 0.0
            out_boundary = xind > (nx - 1.5)
            mask[out_boundary] = 0.0  # Set to zero if near the boundary
            xind[out_boundary] = nx - 1

        if show and plotting_available:
            plt.show()

        # Set xind to -1 if in the inner boundary, nx if in outer boundary
        in_boundary = xind < 0.5
        xind[in_boundary] = -1
        out_boundary = xind > (nx - 1.5)
        xind[out_boundary] = nx

        return xind.reshape(input_shape), zind.reshape(input_shape)

    def metric(self):
        """Return the metric tensor, dx and dz

        Returns
        -------
        dict
            Dictionary containing:
            - **dx, dz**: Grid spacing
            - **gxx, gxz, gzz**: Covariant components
            - **g_xx, g_xz, g_zz**: Contravariant components

        """

        dx = 1.0 / float(self.nx - 1)  # x from 0 to 1
        dz = 2.0 * np.pi / float(self.nz)  # z from 0 to 2pi

        # Get arrays of indices
        xind, zind = np.meshgrid(np.arange(self.nx), np.arange(self.nz), indexing="ij")

        # Calculate the gradient along each coordinate
        dRdx, dZdx = self.getCoordinate(xind, zind, dx=1)
        dRdx /= dx
        dZdx /= dx
        dRdz, dZdz = self.getCoordinate(xind, zind, dz=1)
        dRdz /= dz
        dZdz /= dz

        g_xx = dRdx ** 2 + dZdx ** 2
        g_xz = dRdx * dRdz + dZdx * dZdz
        g_zz = dRdz ** 2 + dZdz ** 2

        # Calculate metric by inverting
        # ( gxx   gxz ) = ( g_xx   g_xz )^-1
        # ( gxz   gzz )   ( g_xz   g_zz )

        determinant = g_xx * g_zz - g_xz ** 2
        gxx = g_zz / determinant
        gzz = g_xx / determinant
        gxz = -g_xz / determinant

        return {
            "dx": dx,
            "dz": dz,  # Grid spacing
            "gxx": gxx,
            "g_xx": g_xx,
            "gxz": gxz,
            "g_xz": g_xz,
            "gzz": gzz,
            "g_zz": g_zz,
        }


def grid_annulus(inner, outer, nx, nz, show=True, return_coords=False):
    """Grid an annular region, given inner and outer boundaries both of
    which are RZline objects

    This is a very simple algorithm which just draws straight lines
    between inner and outer boundaries.

    Parameters
    ----------
    inner, outer : `RZline`
        Inner and outer boundaries of the domain
    nx : int
        The required radial resolution, including boundaries
    nz : int
        The required poloidal resolution
    show : bool, optional
        If True, plot the resulting grid
    return_coords : bool, optional
        If True, return the R, Z coordinates of the grid points,
        instead of a `StructuredPoloidalGrid`

    Returns
    -------
    StructuredPoloidalGrid
        A grid of the region

    """

    assert nx >= 2
    assert nz > 1

    R = zeros((nx, nz))
    Z = zeros((nx, nz))

    # Generate angle values, which should now be equally spaced
    # in distance along inner and outer boundaries
    thetavals = linspace(0, 2 * pi, nz, endpoint=False)

    # Radial coordinate
    xvals = linspace(0, 1.0, nx, endpoint=True)

    innerR = inner.Rvalue(thetavals)
    innerZ = inner.Zvalue(thetavals)

    outerR = outer.Rvalue(thetavals)
    outerZ = outer.Zvalue(thetavals)
    for i, x in enumerate(xvals):
        # Get the R and Z coordinates of this line
        R[i, :] = x * outerR + (1.0 - x) * innerR
        Z[i, :] = x * outerZ + (1.0 - x) * innerZ

    if show and plotting_available:
        plt.plot(inner.R, inner.Z, "-o")
        plt.plot(outer.R, outer.Z, "-o")

        plt.plot(R, Z, "x")

        plt.show()

    if return_coords:
        return R, Z
    return StructuredPoloidalGrid(R, Z)


def grid_elliptic(
    inner,
    outer,
    nx,
    nz,
    show=False,
    tol=1e-10,
    align=True,
    restrict_size=20,
    restrict_factor=2,
    return_coords=False,
):
    """Create a structured grid between inner and outer boundaries using
    elliptic method

    Coordinates x = x(R, Z) and z = z(R,Z) obey an elliptic equation:

    .. math::

        d^2x/dR^2 + d^2x/dZ^2 = 0

        d^2z/dR^2 + d^2z/dZ^2 = 0

    where here x is in in the domain (0, 1) and z in (0, 2pi)

    The above equations are inverted, giving:

    .. math::

        a*R_xx - 2*b*R_xz + c*R_zz = 0

        a*Z_xx - 2*b*Z_xz + c*Z_zz = 0

    where

    .. math::

        a &= R_z^2 + Z_z^2

        b &= R_z*R_x + Z_x*Z_z

        c &= R_x^2 + Z_x^2

    This is a nonlinear system of equations which is solved
    iteratively.

    Parameters
    ----------
    inner, outer : `RZline`
        Inner and outer boundaries of the domain
    nx : int
        The required radial resolution, including boundaries
    nz : int
        The required poloidal resolution
    show : bool, optional
        Display plots of intermediate results
    tol : float, optional
        Controls when iteration stops
    align : bool, optional
        Attempt to align the inner and outer boundaries
    restrict_size : int, optional
        The size (nx or nz) above which the grid is coarsened
    restrict_factor : int, optional
        The factor by which the grid is divided if coarsened
    return_coords : bool, optional
        If True, return the R, Z coordinates of the grid points,
        instead of a `StructuredPoloidalGrid`

    Returns
    -------
    If return_coords is true, returns R,Z as arrays.
    If return_coords is false, returns a `StructuredPoloidalGrid` object

    References
    ----------
    https://www.nada.kth.se/kurser/kth/2D1263/l2.pdf
    https://en.wikipedia.org/wiki/Principles_of_grid_generation

    """

    assert nx >= 2
    assert nz > 1

    # Generate angle values (y coordinate),
    # which should now be equally spaced
    # in distance along inner and outer boundaries
    thetavals = linspace(0, 2 * pi, nz, endpoint=False)

    # Radial coordinate
    xvals = linspace(0, 1.0, nx, endpoint=True)

    if align:
        # Align inner and outer boundaries
        # Easiest way is to roll both boundaries
        # so that index 0 is on the outboard midplane

        ind = np.argmax(inner.R)
        inner = rzline.RZline(np.roll(inner.R, -ind), np.roll(inner.Z, -ind))
        ind = np.argmax(outer.R)
        outer = rzline.RZline(np.roll(outer.R, -ind), np.roll(outer.Z, -ind))

    if (nx > restrict_size) or (nz > restrict_size):
        # Create a coarse grid first to get a starting guess
        # Only restrict the dimensions which exceed restrict_size
        # Note that this might result in multiple levels of resolution

        nx_r = nx
        if nx > restrict_size:
            nx_r = int(nx / restrict_factor)

        nz_r = nz
        if nz > restrict_size:
            nz_r = int(nz / restrict_factor)

        # Create the coarse mesh
        R_r, Z_r = grid_elliptic(
            inner,
            outer,
            nx_r,
            nz_r,
            align=False,
            tol=tol,
            restrict_size=restrict_size,
            restrict_factor=restrict_factor,
            return_coords=True,
        )

        # Note: Lower case x,z are indices
        z_r = linspace(
            0, 2 * pi, nz_r + 1, endpoint=True
        )  # Add on the final point duplicating the first
        x_r = linspace(0, 1.0, nx_r, endpoint=True)

        # Note: Upper case R,Z are real-space locations
        R_r = np.concatenate((R_r, np.reshape(R_r[:, 0], (nx_r, 1))), axis=1)
        Z_r = np.concatenate((Z_r, np.reshape(Z_r[:, 0], (nx_r, 1))), axis=1)

        # Now interpolate
        spl = RectBivariateSpline(x_r, z_r, R_r)
        R = spl(xvals, thetavals, grid=True)
        spl = RectBivariateSpline(x_r, z_r, Z_r)
        Z = spl(xvals, thetavals, grid=True)

        # Make sure that the inner and outer boundaries are on the
        # inner and outer RZline, not interpolated
        R[0, :] = inner.Rvalue(thetavals)
        Z[0, :] = inner.Zvalue(thetavals)

        R[-1, :] = outer.Rvalue(thetavals)
        Z[-1, :] = outer.Zvalue(thetavals)

    else:
        # Interpolate coordinates of inner and outer boundary
        Rinner = inner.Rvalue(thetavals)
        Zinner = inner.Zvalue(thetavals)

        Router = outer.Rvalue(thetavals)
        Zouter = outer.Zvalue(thetavals)

        # Interpolate in x between inner and outer
        # to get starting guess for a grid
        R = zeros((nx, nz))
        Z = zeros((nx, nz))
        for i in range(nx):
            R[i, :] = xvals[i] * Router + (1.0 - xvals[i]) * Rinner
            Z[i, :] = xvals[i] * Zouter + (1.0 - xvals[i]) * Zinner

    dx = xvals[1] - xvals[0]
    dz = thetavals[1] - thetavals[0]

    if show and plotting_available:
        # Markers on original points on inner and outer boundaries
        plt.plot(inner.R, inner.Z, "-o")
        plt.plot(outer.R, outer.Z, "-o")

        # Black lines through inner and outer boundaries
        r, z = inner.position(np.linspace(0, 2 * np.pi, 100))
        plt.plot(r, z, "k")
        r, z = outer.position(np.linspace(0, 2 * np.pi, 100))
        plt.plot(r, z, "k")

        # Red dots to mark the inner and outer boundaries
        plt.plot(R[0, :], Z[0, :], "ro")
        plt.plot(R[-1, :], Z[-1, :], "ro")

    # Start solver loop
    while True:
        # Calculate coefficients, which exclude boundary points
        # Note that the domain is periodic in y so roll arrays

        R_xm = R[:-2, :]  # R(x-1,z)
        R_xp = R[2:, :]  # R(x+1,z)
        R_zm = np.roll(R, 1, axis=1)  # R(x, z-1)
        R_zp = np.roll(R, -1, axis=1)  # R(x, z+1)
        R_xmzm = R_zm[:-2, :]  # R(x-1, z-1)
        R_xpzm = R_zm[2:, :]  # R(x+1, z-1)
        R_xmzp = R_zp[:-2, :]  # R(x-1, z+1)
        R_xpzp = R_zp[2:, :]  # R(x+1, z+1)
        R_zm = R_zm[1:-1, :]  # Now chop off x boundaries
        R_zp = R_zp[1:-1, :]  # This is done to minimise number of rolls

        Z_xm = Z[:-2, :]
        Z_xp = Z[2:, :]
        Z_zm = np.roll(Z, 1, axis=1)
        Z_zp = np.roll(Z, -1, axis=1)
        Z_xmzm = Z_zm[:-2, :]
        Z_xpzm = Z_zm[2:, :]
        Z_xmzp = Z_zp[:-2, :]
        Z_xpzp = Z_zp[2:, :]
        Z_zm = Z_zm[1:-1, :]
        Z_zp = Z_zp[1:-1, :]

        dRdz = (R_zp - R_zm) / (2.0 * dz)
        dRdx = (R_xp - R_xm) / (2.0 * dx)

        dZdz = (Z_zp - Z_zm) / (2.0 * dz)
        dZdx = (Z_xp - Z_xm) / (2.0 * dx)

        a = dRdz ** 2 + dZdz ** 2
        b = dRdz * dRdx + dZdx * dZdz
        c = dRdx ** 2 + dZdx ** 2

        # Now solve a*R_xx - 2*b*R_xz + c*R_zz = 0
        # For now using Jacobi update

        a_dx2 = a / dx ** 2
        b_dxdz = b / (2.0 * dx * dz)
        c_dz2 = c / dz ** 2
        inv_diag = 1.0 / (2 * a / dx ** 2 + 2 * c / dz ** 2)

        Rold = R.copy()
        Zold = Z.copy()

        R[1:-1, :] = (
            a_dx2 * (R_xm + R_xp)
            - b_dxdz * (R_xpzp - R_xmzp - R_xpzm + R_xmzm)
            + c_dz2 * (R_zm + R_zp)
        ) * inv_diag

        Z[1:-1, :] = (
            a_dx2 * (Z_xm + Z_xp)
            - b_dxdz * (Z_xpzp - Z_xmzp - Z_xpzm + Z_xmzm)
            + c_dz2 * (Z_zm + Z_zp)
        ) * inv_diag

        maxchange_sq = np.amax((R - Rold) ** 2 + (Z - Zold) ** 2)

        if maxchange_sq < tol:
            break

    if show and plotting_available:
        plt.plot(R, Z)
        plt.plot(np.transpose(R), np.transpose(Z))
        plt.show()

    if return_coords:
        return R, Z
    return StructuredPoloidalGrid(R, Z)


if __name__ == "__main__":

    import rzline

    # inner = circle(R0=1.5, r=1.0, n=100)
    # outer = circle(R0=1.0, r=2.0, n=100)

    inner = rzline.shaped_line(R0=3.0, a=0.5, elong=1.0, triang=0.0, indent=1.0, n=50)
    outer = rzline.shaped_line(R0=2.8, a=1.5, elong=1.0, triang=0.0, indent=0.2, n=50)
    # outer = shaped_line(R0=3.0, a=1.0, elong=1.0, triang=0.0, indent=1.0, n=50)

    grid = grid_elliptic(inner, outer, 100, 100, show=True)

    # grid.findIndex(2.0, 1.5)
    x, z = grid.findIndex([2.0, 1.9], [1.5, 2.0])
    print(x, z)
