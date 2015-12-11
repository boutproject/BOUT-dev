from __future__ import division
from builtins import object
from past.utils import old_div

import numpy as np

class Grid(object):
    def __init__(self, nx, ny, nz,
                 Lx=0.1, Ly=10., Lz = 1.,
                 name="fci_grid", MXG=2):
        """Initialise a grid object

        Inputs
        ------
        nx  - Number of radial points
        ny  - Number of toroidal points (NOTE: Different to BOUT++ standard)
        nz  - Number of poloidal points

        Lx  - Radial domain size  [m]
        Ly  - Toroidal domain size [m]
        Lz  - Poloidal domain size [m]

        name - Grid name (default "fci_grid")

        MXG - Number of x guard cells
        """

        self.MXG = MXG

        self.name = name

        self.nx = int(nx)
        self.ny = int(ny)
        self.nz = int(nz)

        self.Lx = float(Lx)
        self.Ly = float(Ly)
        self.Lz = float(Lz)

        self.delta_x = old_div(self.Lx,(nx-2.*MXG))
        self.delta_y = old_div(self.Ly,ny)
        self.delta_z = old_div(self.Lz,nz)

        # Coord arrays
        self.xarray = Lx * (np.arange(nx) - MXG + 0.5)/(nx - 2.*MXG)  # 0 and 1 half-way between cells
        self.yarray = np.linspace(0,Ly,ny,endpoint=False)
        self.zarray = np.linspace(0,Lz,nz,endpoint=False)

        self.xcentre = old_div(Lx, 2.)#0.5*max(self.xarray)
        self.zcentre = 0.5*max(self.zarray)

        # How to do this properly?
        def Rmaj(x,z,phi):
            return 1.0

        self.Rmaj = Rmaj

    def __repr__(self):
        return "Grid(nx={nx}, ny={ny}, nz={nz}, Lx={Lx}, Ly={Ly}, Lz={Lz}, name='{name}')".format(
            nx=self.nx, ny=self.ny, nz=self.nz, Lx=self.Lx, Ly=self.Ly, Lz=self.Lz, name=self.name)
