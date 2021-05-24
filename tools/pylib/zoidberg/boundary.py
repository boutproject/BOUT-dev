"""Boundary objects that define an 'outside'

"""

import numpy as np


class NoBoundary(object):
    """No boundary, so no points outside"""

    def outside(self, x, y, z):
        """Returns True if the point is outside the boundary

        Parameters
        ----------
        x, y, z : array_like
            Coordinates of the point(s) to check

        Returns
        -------
        bool
            True if point is outside boundary
        """
        return np.full(x.shape, False, dtype=bool)


class RectangularBoundaryXZ(object):
    def __init__(self, xmin, xmax, zmin, zmax):
        """Create a simple rectangular boundary box

        Parameters
        ----------
        xmin, xmax : float
            The X coordinates of the limits of the box
        zmin, zmax : float
            The Z coordinates of the limits of the box
        """

        self.xmin = xmin
        self.xmax = xmax
        self.zmin = zmin
        self.zmax = zmax

    def outside(self, x, y, z):
        """Returns true if the given point is outside the domain

        Parameters
        ----------
        x, y, z : array_like
            Coordinates of the point(s) to check

        Returns
        -------
        bool
            True if point is outside boundary

        """

        return np.logical_or(
            np.logical_or(x < self.xmin, x > self.xmax),
            np.logical_or(z < self.zmin, z > self.zmax),
        )


class PolygonBoundaryXZ(object):
    def __init__(self, xarr, zarr):
        """Create a polygon boundary defined by

        (xarr[0], zarr[0]) -- (xarr[1], zarr[1]) -- ...

        Parameters
        ----------
        xarr, zarr : array_like
            The X, Z coordinates of the polygon vertices
        """

        # Find a point outside the polygon
        self.zout = np.amax(zarr) + 1.0  # Somewhere outside the polygon
        self.xarr = xarr
        self.zarr = zarr

    def outside(self, x, y, z):
        """Returns true if the given point is outside the domain

        Parameters
        ----------
        x, y, z : array_like
            Coordinates of the point(s) to check

        Returns
        -------
        bool
            True if point is outside boundary

        """
        # Find intersection of a line from z = z_0 to z_1
        #
        #    x = x_0    z = (1-a)*z_0 + a*z_1
        #
        # with a line segment (x_0, z_0) to (x_1, z_1)
        #
        #    x' = (1-b)*x'_0 + b*x'_1   z' = (1-b)*z'_0 + b*z'_1
        #
        # where a and b must be between 0 and 1 for an intersection
        #

        x_0 = x
        z_0 = z
        z_1 = self.zout

        # Create an array of True values, marking all points
        # as outside the domain
        outside = np.full(x.shape, True, dtype=bool)

        # Each intersection with a boundary segment
        # changes True->False->True i.e. XOR

        nsegments = len(self.xarr)
        for i in range(nsegments):  # Loop over all segments

            xp_0 = self.xarr[i]
            zp_0 = self.zarr[i]

            ip = (i + 1) % nsegments  # Next point
            xp_1 = self.xarr[ip]
            zp_1 = self.zarr[ip]

            # First calculate b = (x_0 - x'_0) / (x'_1 - x'_0)
            # but check for case where (x'_1 - x'_0) is small
            if abs(xp_1 - xp_0) < 1e-6:
                continue  # Skip this segment

            b = (x_0 - xp_0) / (xp_1 - xp_0)

            a = (zp_0 - z_0 + b * (zp_1 - zp_0)) / (z_1 - z_0)

            intersect = np.logical_and(
                np.logical_and(b >= 0.0, b <= 1.0), np.logical_and(a >= 0.0, a <= 1.0)
            )

            # intersect is now true for those points which intersect this line segment
            outside = np.logical_xor(outside, intersect)

        return outside
