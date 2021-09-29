"""Routines and classes for representing periodic lines in R-Z
poloidal planes

"""

import numpy as np

from numpy import pi, linspace, sqrt, cos, sin, append, argmin
from scipy.interpolate import splrep, splev, interp1d
from scipy.integrate import cumtrapz

import warnings

try:
    import matplotlib.pyplot as plt

    plotting_available = True
except ImportError:
    warnings.warn("Couldn't import matplotlib, plotting not available.")
    plotting_available = False


class RZline:
    """Represents (R,Z) coordinates of a periodic line

    Attributes
    ----------
    R : array_like
        Major radius [m]
    Z : array_like
        Height [m]
    theta : array_like
        Angle variable [radians]

        `R`, `Z` and `theta` all have the same length

    Parameters
    ----------
    r, z : array_like
        1D arrays of the major radius (`r`) and height (`z`) which are
        of the same length. A periodic domain is assumed, so the last
        point connects to the first.

    anticlockwise : bool, optional
        Ensure that the line goes anticlockwise in the R-Z plane
        (positive theta)

    Note that the last point in (r,z) arrays should not be the same
    as the first point. The (r,z) points are in [0,2pi)

    The input r,z points will be reordered, so that the
    theta angle goes anticlockwise in the R-Z plane

    """

    def __init__(self, r, z, anticlockwise=True):
        r = np.asfarray(r)
        z = np.asfarray(z)

        # Check the sizes of the variables
        n = len(r)
        assert len(z) == n
        assert r.shape == (n,)
        assert z.shape == (n,)

        if anticlockwise:
            # Ensure that the line is going anticlockwise (positive theta)
            mid_ind = np.argmax(r)  # Outboard midplane index
            if z[(mid_ind + 1) % n] < z[mid_ind]:
                # Line going down at outboard midplane. Need to reverse
                r = r[::-1]  # r = np.flip(r)
                z = z[::-1]  # z = np.flip(z)

        self.R = r
        self.Z = z

        # Define an angle variable
        self.theta = linspace(0, 2 * pi, n, endpoint=False)

        # Create a spline representation
        # Note that the last point needs to be passed but is not used
        self._rspl = splrep(append(self.theta, 2 * pi), append(r, r[0]), per=True)
        self._zspl = splrep(append(self.theta, 2 * pi), append(z, z[0]), per=True)

    def Rvalue(self, theta=None, deriv=0):
        """Calculate the value of R at given theta locations

        Parameters
        ----------
        theta : array_like, optional
            Theta locations to find R at. If None (default), use the
            values of theta stored in the instance
        deriv : int, optional
            The order of derivative to compute (default is just the R value)

        Returns
        -------
        ndarray
            Value of R at each input theta point
        """
        if theta is None:
            theta = self.theta
        else:
            theta = np.remainder(theta, 2 * np.pi)

        return splev(theta, self._rspl, der=deriv)

    def Zvalue(self, theta=None, deriv=0):
        """Calculate the value of Z at given theta locations

        Parameters
        ----------
        theta : array_like, optional
            Theta locations to find Z at. If None (default), use the
            values of theta stored in the instance
        deriv : int, optional
            The order of derivative to compute (default is just the Z value)

        Returns
        -------
        ndarray
            Value of Z at each input theta point
        """
        if theta is None:
            theta = self.theta
        else:
            theta = np.remainder(theta, 2 * np.pi)
        return splev(theta, self._zspl, der=deriv)

    def position(self, theta=None):
        """Calculate the value of both R, Z at given theta locations

        Parameters
        ----------
        theta : array_like, optional
            Theta locations to find R, Z at. If None (default), use the
            values of theta stored in the instance

        Returns
        -------
        R, Z : (ndarray, ndarray)
            Value of R, Z at each input theta point
        """
        return self.Rvalue(theta=theta), self.Zvalue(theta=theta)

    def positionPolygon(self, theta=None):
        """Calculates (R,Z) position at given theta angle by joining points
        by straight lines rather than a spline. This avoids the
        overshoots which can occur with splines.

        Parameters
        ----------
        theta : array_like, optional
            Theta locations to find R, Z at. If None (default), use the
            values of theta stored in the instance

        Returns
        -------
        R, Z : (ndarray, ndarray)
            Value of R, Z at each input theta point
        """
        if theta is None:
            return self.R, self.Z
        n = len(self.R)
        theta = np.remainder(theta, 2.0 * pi)
        dtheta = 2.0 * np.pi / n
        ind = np.trunc(theta / dtheta)
        rem = np.remainder(theta, dtheta)
        indp = (ind + 1) % n
        return (rem * self.R[indp] + (1.0 - rem) * self.R[ind]), (
            rem * self.Z[indp] + (1.0 - rem) * self.Z[ind]
        )

    def distance(self, sample=20):
        """Integrates the distance along the line.

        Parameters
        ----------
        sample : int, optional
            Number of samples to take per point

        Returns
        -------
        An array one longer than theta. The first element is zero,
        and the last element is the total distance around the loop

        """

        sample = int(sample)
        assert sample >= 1

        thetavals = np.linspace(
            0.0, 2.0 * np.pi, sample * len(self.theta) + 1, endpoint=True
        )

        # variation of length with angle dl/dtheta
        dldtheta = sqrt(
            self.Rvalue(thetavals, deriv=1) ** 2 + self.Zvalue(thetavals, deriv=1) ** 2
        )

        # Integrate cumulatively, then take only the values at the grid points (including end)
        return cumtrapz(dldtheta, thetavals, initial=0.0)[::sample]

    def equallySpaced(self, n=None):
        """Returns a new RZline which has a theta uniform in distance along
        the line

        Parameters
        ----------
        n : int, optional
            Number of points. Default is the same as the current line

        Returns
        -------
        `RZline`
            A new `RZline` based on this instance, but with uniform theta-spacing
        """
        if n is None:
            n = len(self.theta)

        # Distance along the line
        dist = self.distance()

        # Positions where points are desired
        positions = linspace(dist[0], dist[-1], n, endpoint=False)

        # Find which theta value these correspond to
        thetavals = interp1d(
            dist, append(self.theta, 2.0 * pi), copy=False, assume_sorted=True
        )
        new_theta = thetavals(positions)

        return RZline(self.Rvalue(new_theta), self.Zvalue(new_theta))

    def closestPoint(self, R, Z, niter=3, subdivide=20):
        """Find the closest point on the curve to the given (R,Z) point

        Parameters
        ----------
        R, Z : float
            The input R, Z point
        niter : int, optional
            How many iterations to use

        Returns
        -------
        float
            The value of theta (angle)

        """

        # First find the closest control point
        ind = argmin((self.R - R) ** 2 + (self.Z - Z) ** 2)
        theta0 = self.theta[ind]
        dtheta = self.theta[1] - self.theta[0]

        # Iteratively refine and find new minimum
        for i in range(niter):
            # Create a new set of points between point (ind +/- 1)
            # By using dtheta, wrapping around [0,2pi] is handled
            thetas = np.linspace(
                theta0 - dtheta, theta0 + dtheta, subdivide, endpoint=False
            )
            Rpos, Zpos = self.positionPolygon(thetas)

            ind = argmin((Rpos - R) ** 2 + (Zpos - Z) ** 2)
            theta0 = thetas[ind]
            dtheta = thetas[1] - thetas[0]

        return np.remainder(theta0, 2 * np.pi)

    def plot(self, axis=None, show=True):
        """Plot the RZline, either on the given axis or a new figure

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

        theta = np.linspace(0, 2 * np.pi, 100, endpoint=True)
        axis.plot(self.Rvalue(theta), self.Zvalue(theta), "k-")
        axis.plot(self.R, self.Z, "ro")

        if show:
            plt.show()

        return axis


def circle(R0=1.0, r=0.5, n=20):
    """Creates a pair of RZline objects, for inner and outer boundaries

    Parameters
    ----------
    R0 : float, optional
        Centre point of the circle
    r : float, optional
        Radius of the circle
    n : int, optional
        Number of points to use in the boundary

    Returns
    -------
    `RZline`
        A circular `RZline`

    """
    # Define an angle coordinate
    theta = linspace(0, 2 * pi, n, endpoint=False)

    return RZline(R0 + r * cos(theta), r * sin(theta))


def shaped_line(R0=3.0, a=1.0, elong=0.0, triang=0.0, indent=0.0, n=20):
    """Parametrisation of plasma shape from J. Manickam, Nucl. Fusion 24
    595 (1984)

    Parameters
    ----------
    R0 : float, optional
        Major radius
    a : float, optional
        Minor radius
    elong : float, optional
        Elongation, 0 for a circle
    triang : float, optional
        Triangularity, 0 for a circle
    indent : float, optional
        Indentation, 0 for a circle

    Returns
    -------
    `RZline`
        An `RZline` matching the given parameterisation

    """
    theta = linspace(0, 2 * pi, n, endpoint=False)
    return RZline(
        R0 - indent + (a + indent * cos(theta)) * cos(theta + triang * sin(theta)),
        (1.0 + elong) * a * sin(theta),
    )


def line_from_points_poly(rarray, zarray, show=False):
    """Find a periodic line which goes through the given (r,z) points

    This function starts with a triangle, then adds points
    one by one, inserting into the polygon along the nearest
    edge

    Parameters
    ----------
    rarray, zarray : array_like
        R, Z coordinates. These arrays should be the same length

    Returns
    -------
    `RZline`
        An `RZline` object representing a periodic line

    """

    rarray = np.asfarray(rarray)
    zarray = np.asfarray(zarray)

    assert rarray.size >= 3
    assert rarray.shape == zarray.shape

    npoints = rarray.size

    rvals = rarray.copy()
    zvals = zarray.copy()

    # Take the first three points to make a triangle

    if show and plotting_available:
        plt.figure()
        plt.plot(rarray, zarray, "x")
        plt.plot(
            np.append(rvals[:3], rvals[0]), np.append(zvals[:3], zvals[0])
        )  # Starting triangle

    for i in range(3, npoints):
        line = RZline(rvals[:i], zvals[:i])

        angle = np.linspace(0, 2 * pi, 100)
        r, z = line.position(angle)

        # Next point to add
        # plt.plot(rarray[i], zarray[i], 'o')

        # Find the closest point on the line
        theta = line.closestPoint(rarray[i], zarray[i])

        rl, zl = line.position(theta)

        ind = np.floor(float(i) * theta / (2.0 * np.pi))

        # Insert after this index

        if ind != i - 1:
            # If not the last point, then need to shift other points along
            rvals[ind + 2 : i + 1] = rvals[ind + 1 : i]
            zvals[ind + 2 : i + 1] = zvals[ind + 1 : i]
        rvals[ind + 1] = rarray[i]
        zvals[ind + 1] = zarray[i]

        if show and plotting_available:
            plt.plot([rarray[i], rl], [zarray[i], zl])
            plt.plot(
                np.append(rvals[: (i + 1)], rvals[0]),
                np.append(zvals[: (i + 1)], zvals[0]),
            )  # New line

    if show and plotting_available:
        plt.show()
    return RZline(rvals, zvals)


def line_from_points(rarray, zarray, show=False):
    """Find a periodic line which goes through the given (r,z) points

    This function starts at a point, and finds the nearest neighbour
    which is not already in the line

    Parameters
    ----------
    rarray, zarray : array_like
        R, Z coordinates. These arrays should be the same length

    Returns
    -------
    `RZline`
        An `RZline` object representing a periodic line

    """

    # Make sure we have Numpy arrays
    rarray = np.asfarray(rarray)
    zarray = np.asfarray(zarray)

    assert rarray.size == zarray.size

    # We can get different answers depending on which point
    # we start the line on.
    # Therefore start the line from every point in turn,
    # and keep the line with the shortest total distance

    best_line = None  # The best line found so far
    best_dist = 0.0  # Distance around best line

    for start_ind in range(rarray.size):

        # Create an array of remaining points
        # Make copies since we edit the array later
        rarr = np.roll(rarray, start_ind).copy()
        zarr = np.roll(zarray, start_ind).copy()

        # Create new lists for the result
        rvals = [rarr[0]]
        zvals = [zarr[0]]

        rarr = rarr[1:]
        zarr = zarr[1:]

        while rarr.size > 1:
            # Find the index in array closest to last point
            ind = np.argmin((rvals[-1] - rarr) ** 2 + (zvals[-1] - zarr) ** 2)

            rvals.append(rarr[ind])
            zvals.append(zarr[ind])
            # Shift arrays
            rarr[ind:-1] = rarr[(ind + 1) :]
            zarr[ind:-1] = zarr[(ind + 1) :]
            # Chop off last point
            rarr = rarr[:-1]
            zarr = zarr[:-1]

        # One left, add to the end
        rvals.append(rarr[0])
        zvals.append(zarr[0])

        new_line = RZline(rvals, zvals)
        new_dist = new_line.distance()[-1]  # Total distance

        if (best_line is None) or (new_dist < best_dist):
            # Either if we haven't got a line, or found
            # a better line
            best_line = new_line
            best_dist = new_dist

    return best_line


if __name__ == "__main__":

    import field
    import fieldtracer
    import poloidal_grid

    #############################################################################
    # Define the magnetic field

    # Length in y after which the coils return to their starting (R,Z) locations
    yperiod = 10.0

    magnetic_field = field.StraightStellarator(I_coil=0.3, radius=1.0, yperiod=yperiod)

    #############################################################################
    # Create the inner flux surface, starting at a point at phi=0
    # To do this we need to define the y locations of the poloidal points
    # where we will construct grids

    start_r = 0.2
    start_z = 0.0

    nslices = 8  # Number of poloidal slices
    ycoords = np.linspace(0, yperiod, nslices)
    npoints = 20  # Points per poloidal slice

    # Create a field line tracer
    tracer = fieldtracer.FieldTracer(magnetic_field)

    # Extend the y coordinates so the tracer loops npoints times around yperiod
    ycoords_all = ycoords
    for i in range(1, npoints):
        ycoords_all = np.append(ycoords_all, ycoords + i * yperiod)

    coord = tracer.follow_field_lines(start_r, start_z, ycoords_all, rtol=1e-12)

    inner_lines = []
    for i in range(nslices):
        r = coord[i::nslices, 0]
        z = coord[i::nslices, 1]
        line = line_from_points(r, z)
        # Re-map the points so they're approximately uniform in distance along the surface
        # Note that this results in some motion of the line
        line = line.equallySpaced()
        inner_lines.append(line)

    # Now have a list of y coordinates (ycoords) and inner lines (inner_lines)

    #############################################################################
    # Generate a fixed circle for the outer boundary

    outer_line = circle(R0=0.0, r=0.8)

    #############################################################################
    # Now have inner and outer boundaries for each poloidal slice
    # Generate a grid on each poloidal slice using the elliptic grid generator

    nx = 20
    ny = 20

    pol_slices = [
        poloidal_grid.grid_elliptic(inner_line, outer_line, nx, ny, show=True)
        for inner_line in inner_lines
    ]
