#!/usr/bin/env python3
from __future__ import print_function
from builtins import range

import numpy as np
from numpy import cos, sin, pi

from tvtk.api import tvtk
#from enthought.mayavi.scripts import mayavi2

def aligned_points(grid, nz=1, period=1.0, maxshift=0.4):
    try:
        nx = grid["nx"]#[0]
        ny = grid["ny"]#[0]
        zshift = grid["zShift"]
        Rxy = grid["Rxy"]
        Zxy = grid["Zxy"]
    except:
        print("Missing required data")
        return None
    
    
    dz = 2.*pi / (period * (nz-1))
    phi0 = np.linspace(0,2.*pi / period, nz)
    
    
    # Need to insert additional points in Y so mesh looks smooth
    #for y in range(1,ny):
    #    ms = np.max(np.abs(zshift[:,y] - zshift[:,y-1]))
    #    if(

    # Create array of points, structured

    points = np.zeros([nx*ny*nz, 3])
    
    
    start = 0
    for y in range(ny):
        
       
        end = start + nx*nz
        
        phi = zshift[:,y] + phi0[:,None]
        r = Rxy[:,y] + (np.zeros([nz]))[:,None]
    
        xz_points = points[start:end]
    
        
        xz_points[:,0] = (r*cos(phi)).ravel()  # X
        xz_points[:,1] = (r*sin(phi)).ravel()  # Y
        xz_points[:,2] = (Zxy[:,y]+(np.zeros([nz]))[:,None]).ravel()   # Z
    
        
        start = end
    
    return points

def create_grid(grid, data, period=1):
    
    s = np.shape(data)
    
    nx = grid["nx"]#[0]
    ny = grid["ny"]#[0]
    nz = s[2]
    
    print("data: %d,%d,%d   grid: %d,%d\n" % (s[0],s[1],s[2], nx,ny))
    
    dims = (nx, nz, ny)
    sgrid = tvtk.StructuredGrid(dimensions=dims)
    pts = aligned_points(grid, nz, period)
    print(np.shape(pts))
    sgrid.points = pts
    
    scalar = np.zeros([nx*ny*nz])
    start = 0
    for y in range(ny):
        end = start + nx*nz
        
        #scalar[start:end] = (data[:,y,:]).transpose().ravel()
        scalar[start:end] = (data[:,y,:]).ravel()

        print(y, " = " , np.max(scalar[start:end]))
        start = end
    
    sgrid.point_data.scalars = np.ravel(scalar.copy())
    sgrid.point_data.scalars.name = "data"
    
    return sgrid

#@mayavi2.standalone
def view3d(sgrid):
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi.modules.api import Outline, GridPlane
    from mayavi.api import Engine
    from mayavi.core.ui.engine_view import EngineView
    e=Engine()
    e.start()
    s = e.new_scene()
     # Do this if you need to see the MayaVi tree view UI.
    ev = EngineView(engine=e)
    ui = ev.edit_traits()

#    mayavi.new_scene()
    src = VTKDataSource(data=sgrid)
    e.add_source(src)
    e.add_module(Outline())
    g = GridPlane()
    g.grid_plane.axis = 'x'
    e.add_module(g)

if __name__ == '__main__':
    from boutdata.collect import collect
    from boututils.file_import import file_import
    
    #path = "/media/449db594-b2fe-4171-9e79-2d9b76ac69b6/runs/data_33/"
    path="../data"

    g = file_import("../bout.grd.nc")
    #g = file_import("../cbm18_8_y064_x516_090309.nc")
    #g = file_import("/home/ben/run4/reduced_y064_x256.nc")
    
    data = collect("P", tind=10, path=path)
    data = data[0,:,:,:]
    s = np.shape(data)
    nz = s[2]
    
    #bkgd = collect("P0", path=path)
    #for z in range(nz):
     #   data[:,:,z] += bkgd

    # Create a structured grid
    sgrid = create_grid(g, data, 1)
    

    w = tvtk.XMLStructuredGridWriter(input=sgrid, file_name='sgrid.vts')
    w.write()
    
    # View the structured grid
    view3d(sgrid)
