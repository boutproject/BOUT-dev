import zoidberg as zb
import numpy as np
import sys
import boutconfig as bc


def rotating_ellipse(
    nx=68,
    ny=16,
    nz=128,
    npoints=421,
    xcentre=5.5,
    I_coil=0.01,
    curvilinear=True,
    rectangular=False,
    fname="rotating-ellipse.fci.nc",
    a=0.4,
    Btor=2.5,
):
    yperiod = 2 * np.pi / 5.0
    field = zb.field.RotatingEllipse(
        xcentre=xcentre,
        I_coil=I_coil,
        radius=2 * a,
        yperiod=yperiod,
        Btor=Btor,
    )
    # Define the y locations
    ycoords = np.linspace(0.0, yperiod, ny, endpoint=False)

    if rectangular:
        print("Making rectangular poloidal grid")
        poloidal_grid = zb.poloidal_grid.RectangularPoloidalGrid(
            nx, nz, 1.0, 1.0, Rcentre=xcentre
        )
    elif curvilinear:
        print("Making curvilinear poloidal grid")
        inner = zb.rzline.shaped_line(
            R0=xcentre, a=a / 2.0, elong=0, triang=0.0, indent=0, n=npoints
        )
        outer = zb.rzline.shaped_line(
            R0=xcentre, a=a, elong=0, triang=0.0, indent=0, n=npoints
        )

        print("creating grid...")
        poloidal_grid = zb.poloidal_grid.grid_elliptic(inner, outer, nx, nz)

    # Create the 3D grid by putting together 2D poloidal grids
    grid = zb.grid.Grid(poloidal_grid, ycoords, yperiod, yperiodic=True)
    maps = zb.make_maps(grid, field, quiet=True)
    zb.write_maps(grid, field, maps, str(fname), metric2d=bc.isMetric2D())


if __name__ == "__main__":
    rotating_ellipse(fname=sys.argv[1])
