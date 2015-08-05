# -*- coding: utf-8 -*-
"""
@author: pcn500

vectory.py Draws vectors that are specified.
"""
#Import the relevant libraries
from visboutit import visual, circle
import os
import numpy as np
from math import sin,cos,pi,fabs
#from boutdata import collect
from scipy import interpolate


#Returns a list of unit vectors that can be dotted with the V_parallel to form the vectors
   
#Clocckwise model
def find_y_unit(ny):
#    y_unit = np.empty((ny))
    y_unit = []
    y_unit_x = np.empty((ny))
    y_unit_y = np.empty((ny))
#Find initial y_unit
    delta_phi = 2*pi*1/(ny-1)
    y_unit_initial =  np.dot((1./(delta_phi)) , [-cos(0) + cos(delta_phi) , sin(0) - sin(delta_phi)])
#    
#    y_unit_x[0] = np.dot((1./(delta_phi)) , [cos(0)- cos(delta_phi),sin(0)-sin(delta_phi)])[0]
#    y_unit_y[0] = np.dot((1./(delta_phi)) , [cos(0)- cos(delta_phi),sin(0)-sin(delta_phi)])[1]
#    y_unit_initial_x = y_unit_x[0]
#    y_unit_initial_y = y_unit_y[0]

#Rotate initial points to find the other remaining values    
    for j in range(ny):
        #Calculate the angle transversed
        phi = 2*pi*j/(ny-1)
        #Update the rotation matrix
        r_mtrx = [[cos(phi) , -sin(phi)] , [sin(phi) , cos(phi) ]]
        y_unit.append(np.dot(r_mtrx, y_unit_initial))
#        y_unit_x[j] = np.dot(r_mtrx , y_unit_initial_x)[0]
#        y_unit_Y[j] = np.dot(r_mtrx , y_unit_initial_y)[1]
        
    
    return y_unit






def torus(name, time, step = 0.5, path = None, R = None, r = None , dr = None , Bt = None, q = None , isttok = False, skip = 1):
    """
    Torus function: Convert to Torodial Coordinate system
    
    Inputs
        name: variable name to import (str)
        time: end time slice to convert to (int), if -1 entire dataset is imported
        step: The gap between the y slices for the y interpolation (float < 1) if 0 then raw data displayed.
        path: File path to data, if blank then the current working directory is used to look for data
        
        Torus settings, if left blank defaults are used all are floats.
        R: Major radius
        r: Minor radius
        dr: Radial width of domain
        Bt: Toroidal magnetic field
        q: Safety factor
        
        isttok: Use the ISTTOK specifications, (Boolean)
        
    Outputs
        Timestamped vtk files of the specified variable.
        A text file with the maximum and minimum files.   
    """    

    
    #Get the working dir
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    # Find the max_t, initial value and shape of variable
    max_t, var_0, nx, ny, nz = visual.dim_all(name)
    
    ny_work = ny
   
   #Interpolate setup
    #Get number of new points
    if step != 0:
        ny_work = len(np.arange(0,(ny-1),step))

    #Define empty array with increased number of y values
    var_new = np.empty((nx,ny_work,nz),dtype=float)
    
    if isttok == True: # Use ISTTOK specifications?
        R,r,dr,Bt,q = 0.46, 0.085, 0.02, 0.5, 5
    # Import coordinates from function in circle library
    r,z = circle.coordinates(nx,ny_work,nz,R,r,dr,Bt,q) # toroidal coordinates
    pts = visual.torus(r, z, nx, ny_work, nz)  # vtk grid points
    
    # Set time
    t = time
    # Set t to max_t, if user desires entire dataset
    if t == -1:
        t = max_t
    #If t is greater than number of slices then set to max_t
    if t >= max_t:
        t = max_t
    
    # Make dir for storing vtk files and image files
    #Check if dir for storing vtk files exsists if not make that dir
    if not os.path.exists("batch"):
        os.makedirs("batch")
    #Make dir for storing images
    if not os.path.exists("images"):
        os.makedirs("images")
    
#    #Initial max and min values
#    var_x, var_y, var_z = var_0[:,0,0], var_0[0,:,0], var_0[0,0,:]
#    max_x, max_y, max_z = np.amax(var_x), np.amax(var_y), np.amax(var_z)
#    min_x, min_y, min_z = np.amin(var_x), np.amin(var_y),np.amin(var_z)
#    max = np.amax((max_x, max_y, max_z))
#    min = np.amin((min_x, min_y, min_z))
    
    #Find the y_unit vectors
    y_unit = find_y_unit(ny_work)
    
    q = 0 
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    while q <= t:
        var = visual.var3d(name,q) # collect variable
#        var = np.zeros((nx,ny_work,nz))
#        for i in range(nx):
#            for j in range(ny):
#                for k in range(nz):
#                    var[i,j,k] = 1
        
#        #Find the min and max values
#        var_x, var_y, var_z = var[:,0,0], var[0,:,0], var[0,0,:]
#        max_x, max_y, max_z = np.amax(var_x), np.amax(var_y), np.amax(var_z)
#        min_x, min_y, min_z = np.amin(var_x), np.amin(var_y), np.amin(var_z)
#        max_i = np.amax(np.array((max_x,max_y,max_z)))
#        min_i = np.amin(np.array((min_x,min_y,min_z)))
#        if max_i >= max:
#            max = max_i
#        if min_i <= min:
#            min = min_i
#        
        if ny != ny_work: #Interpolate section for all y values
            for i in range(nx):
                for k in range(nz):
                    var_y = var[i,:,k]
                    y = np.arange(0,ny) # Arrange y values in array
                    f = interpolate.interp1d(y,var_y) # Interpolate the data
                    y_new = np.arange(0,(ny-1),step) # Arrange new y values 
                    var_y_new = f(y_new) # interpolate y values
                    var_new[i,:,k] = var_y_new # Store values in new variable
        else:
            var_new = var
            
        # Create variables for vectors
        vectorx = np.zeros((nx,ny_work,nz),dtype= float)
        vectory = np.zeros((nx,ny_work,nz),dtype= float)
        vectorz = np.zeros((nx,ny_work,nz),dtype= float)
        
#        sign = np.array( ( 0 , 0 ) , dtype = float)
        
        # Assign values
        for i in range(nx):
            for j in range(ny_work):
                for k in range(nz):
#                    for l in range(2): # Find the sign values
#                        if np.abs(y_unit[j])[l] == y_unit[j][l]:
#                            sign[l] = 1
#                        else:
#                            sign[l] = -1
                    vectorx[i,j,k] = var_new[i,j,k] * y_unit[j][0]
                    vectory[i,j,k] = var_new[i,j,k] * y_unit[j][1]
                    vectorz[i,j,k] = var_new[i,j,k] * 0
        
        #Convert coordinate section        
        vrbl = visual.vtk_var(var_new,nx,ny_work,nz) # vtk variable
        vtk_path = visual.write_vtk_vector2(name,pts,vrbl,vectorx,vectory,vectorz,q)
        print "At t = %d, %d Steps remaining" % (q , ((t- q)/skip)) # Progress indicator
        q += skip
    
    #Write the Max and min values to file
#    mm_array = np.array((max,min))
#    np.savetxt('max_min_' + name + '.txt',mm_array)
    return




'''
import data
make vectorx,vectory,vectorz based on import
convert coordinate system and create mesh

how to convert vectors?


write vectors.

'''
