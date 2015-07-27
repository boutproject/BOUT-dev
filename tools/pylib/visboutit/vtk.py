#!/usr/bin/python

'''
vtk.py
======
This python file contains functions required to convert to cylindercal, torodial and ELM systems.
'''

#Import all the libraries
import visual
import sys
import os
import circle
import numpy as np
from scipy import interpolate

'''
ELM: convert to EML coordinate system,
Inputs: variable name and end time
========
Need to add zShift interpolation percent
'''
def elm(name,time,gridfile,zShift_interpolation_percent):
    
    #Get the working dir
    work_dir = os.getcwd()
    grid_file = str(raw_input('Enter the full filename of the gridfile: '))
    grid_dir = work_dir + '/' + grid_file

    # Get the dimensions and initial value of the variable data.    
    max_t, var_0, nx, ny, nz = visual.dim_all(name)
  
    #ELM
    #Import r,z and zshift from grid file
    r = visual.nc_var(grid_dir,'Rxy') # ELM Coordinates
    z = visual.nc_var(grid_dir,'Zxy') # ELM Coordinates
    zshift = visual.nc_var(grid_dir,'zShift') #ELM Coordinates
    
    
    max_z,min_z = visual.z_shift_mm(nx,ny,zshift) #Find the max and minimum zshift values
    z_tol = min_z + (max_z * 0.25) # Set tolerance value ##CHANGE TO IMPORT PERCENT !!##
    
    
    #Interpolate the Rxy and Zxy values
    r2,ny2 = visual.zshift_interp2d(nx,ny,zshift,z_tol,r)
    z2,ny2 = visual.zshift_interp2d(nx,ny,zshift,z_tol,z)
    zshift2,ny2 = visual.zshift_interp2d(nx,ny,zshift,z_tol,zshift)
    
    #Create points for ELM 
    pts2 = visual.elm(r2,z2,zshift2,nx,ny2,nz)
    
    ## Set time
    t = time
    # Set t to max_t, if user desires entire dataset
    if t == -1:
    	t = max_t
    #If t is greater than number of slices then set to max_t
    if t >= max_t:
            t = max_t
    
    ### Make dir for storing vtk files and image files
    #Check if dir for storing vtk files exsists if not make that dir
    if not os.path.exists("batch"):
    	os.makedirs("batch")
    
    #Make dir for storing images
    if not os.path.exists("images"):
    	os.makedirs("images")
    
    #Initial max and min values
    var_x, var_y, var_z = var_0[:,0,0], var_0[0,:,0], var_0[0,0,:]
    max_x, max_y, max_z = np.amax(var_x), np.amax(var_y), np.amax(var_z)
    min_x, min_y, min_z = np.amin(var_x), np.amin(var_y),np.amin(var_z)
    max = np.amax((max_x, max_y, max_z))
    min = np.amin((min_x, min_y, min_z))
    
    #If the variable is pressure, import and add the P0 variable to this
    if name == 'P':
    	p0 = visual.collect("P0") #Import P0
    	add_p = True
    
     
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    for q in range(t):
    	var = visual.var3d(name,q) # collect variable
    
    	#Find the min and max values
    	var_x, var_y, var_z = var[:,0,0], var[0,:,0], var[0,0,:]
    	max_x, max_y, max_z = np.amax(var_x), np.amax(var_y), np.amax(var_z)
    	min_x, min_y, min_z = np.amin(var_x), np.amin(var_y), np.amin(var_z)
    	max_i = np.amax(np.array((max_x,max_y,max_z)))
    	min_i = np.amin(np.array((min_x,min_y,min_z)))
    	if max_i >= max:
    		max = max_i
    	if min_i <= min:
    		min = min_i
    
    	#Adding p0
    	if add_p == True:#Adding P0
                    for i in range (nx): # X value
                            for j in range(ny): # Y Value
                                    var[i,j,:] = var[i,j,:] + p0[i,j]
    
    	
    	var2,ny3 = visual.zshift_interp3d(nx,ny,nz,zshift,z_tol,var)
    
    	'''
    	#Interpolate section for all y values
    	for i in range(nx):
    		for k in range(nz):
    			var_y = var[i,:,k]
    	                y = np.arange(0,ny) # Arrange y values in array
            	        f = interpolate.interp1d(y,var_y) # Interpolate the data
    	                y_new = np.arange(0,(ny-1),step) # Arrange new y values 
            	        var_y_new = f(y_new) # interpolate y values
    	                var_new[i,:,k] = var_y_new # Store values in new variable
    	'''
    	#Convert coordinate section
    	vrbl = visual.vtk_var(var2,nx,ny2,nz) #vtk variable
    #	vrbl = visual.vtk_var(var,nx,ny,nz) # vtk variable
    	vtk_path = visual.write_vtk(name,pts2,vrbl,q,nx,ny2,nz) # write vtk file
    	print "At t = %d, %d Steps remaining" % (q,t-(q+1)) # Progress indicator
    
    #Write the Max and min values to file
    mm_array = np.array((max,min))
    np.savetxt('max_min_' + name + '.txt',mm_array)
    return 

'''
Torus
=========
Torus coordinate system
'''

def torus(name,time,gridfile):
    grid_file = str(raw_input('Enter the full filename of the gridfile: '))
    grid_dir = work_dir + '/' + grid_file
    
    #Get the working dir
    work_dir = os.getcwd()
    ''' 
    not needed
    #Set the working dir to the Parent of the VisBOUTIt folder
    os.chdir(work_dir +'/../')
    #Update the working dir
    work_dir = os.getcwd()
    '''
    # Find the max_t, initial value and shape of variable
    max_t, var_0, nx, ny, nz = visual.dim_all(name)

    #Interpolate setup
    step = 0.5 #interval for new points
    #Get number of new points
    num_new = len(np.arange(0,(ny-1),step))
    #Define empty array with increased number of y values
    var_new = np.empty((nx,num_new,nz),dtype=float)
    
    # Import coordinates from function in circle library
    r,z = circle.coordinates(nx,num_new,nz,R=0.46,r=0.085,dr=0.02,Bt=0.5,q=5) # toroidal coordinates
    pts = visual.torus(r, z, nx, num_new, nz)  # vtk grid points
    
    ## Set time
    t = time
    # Set t to max_t, if user desires entire dataset
    if t == -1:
        t = max_t
    #If t is greater than number of slices then set to max_t
    if t >= max_t:
        t = max_t
    
    ### Make dir for storing vtk files and image files
    #Check if dir for storing vtk files exsists if not make that dir
    if not os.path.exists("batch"):
        os.makedirs("batch")
    #Make dir for storing images
    if not os.path.exists("images"):
        os.makedirs("images")
    
    #Initial max and min values
    var_x, var_y, var_z = var_0[:,0,0], var_0[0,:,0], var_0[0,0,:]
    max_x, max_y, max_z = np.amax(var_x), np.amax(var_y), np.amax(var_z)
    min_x, min_y, min_z = np.amin(var_x), np.amin(var_y),np.amin(var_z)
    max = np.amax((max_x, max_y, max_z))
    min = np.amin((min_x, min_y, min_z))
    
    q = 0 
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    for q in range(t):
    	var = visual.var3d(name,q) # collect variable
    
    	#Find the min and max values
    	var_x, var_y, var_z = var[:,0,0], var[0,:,0], var[0,0,:]
    	max_x, max_y, max_z = np.amax(var_x), np.amax(var_y), np.amax(var_z)
    	min_x, min_y, min_z = np.amin(var_x), np.amin(var_y), np.amin(var_z)
    	max_i = np.amax(np.array((max_x,max_y,max_z)))
    	min_i = np.amin(np.array((min_x,min_y,min_z)))
    	if max_i >= max:
    		max = max_i
    	if min_i <= min:
    		min = min_i
    
    	#Interpolate section for all y values
    	for i in range(nx):
    		for k in range(nz):
    			var_y = var[i,:,k]
    	                y = np.arange(0,ny) # Arrange y values in array
            	        f = interpolate.interp1d(y,var_y) # Interpolate the data
    	                y_new = np.arange(0,(ny-1),step) # Arrange new y values 
            	        var_y_new = f(y_new) # interpolate y values
    	                var_new[i,:,k] = var_y_new # Store values in new variable
    
    	#Convert coordinate section		
    	vrbl = visual.vtk_var(var_new,nx,num_new,nz) # vtk variable
    	vtk_path = visual.write_vtk(name,pts,vrbl,q,nx,num_new,nz) # write vtk file
    	print "At t = %d, %d Steps remaining" % (q,t-(q+1)) # Progress indicator
    
    #Write the Max and min values to file
    mm_array = np.array((max,min))
    np.savetxt('max_min_' + name + '.txt',mm_array)


'''
Cylinder
=======
This function contains the cylinder format
'''
def cylinder(name, time, gridfile):
    #Get the working dir
    work_dir = os.getcwd()
    grid_file = str(raw_input('Enter the full filename of the gridfile: '))
    grid_path = work_dir + '/' + grid_file

    #Find dimensions of data
    max_t, var_0, nx, ny, nz = visual.dim_all(name)

    # Import coordinates from gridfile
    r = visual.nc_var(grid_path,'Rxy') # toroidal coordinates
    z = visual.nc_var(grid_path,'Zxy') # toroidal coordinates
    pts = visual.cylinder(r, z, nx, ny, nz)  # vtk grid points
    
    #Interpolate setup
    step = 0.5 #interval for new points
    #Get number of new points
    num_new = len(np.arange(0,(ny-1),step))
    #Define empty array with increased number of y values
    var_new = np.empty((nx,num_new,nz),dtype=float)
    
    #Set time from input
    t = time
    
    # Set t to max_t, if user desires entire dataset
    if t == -1:
        t = max_t
    #If t is outside time dim of data set to max_t
    if t >= max_t:
        t = max_t
    
    ### Make dir for storing vtk files and image files
    #Check if dir for storing vtk files exsists if not make that dir
    if not os.path.exists("batch"):
        os.makedirs("batch")
    #Make dir for storing images
    if not os.path.exists("images"):
        os.makedirs("images")
    
    
    #Initial max and min values
    var_x, var_y, var_z = var_0[:,0,0], var_0[0,:,0], var_0[0,0,:]
    max_x, max_y, max_z = np.amax(var_x), np.amax(var_y), np.amax(var_z)
    min_x, min_y, min_z = np.amin(var_x), np.amin(var_y),np.amin(var_z)
    max = np.amax((max_x, max_y, max_z))
    min = np.amin((min_x, min_y, min_z))
    
    q = 0
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    for q in range(t):
            var = visual.var3d(name,q) # collect variable
    
            #Find the min and max values
            var_x, var_y, var_z = var[:,0,0], var[0,:,0], var[0,0,:]
            max_x, max_y, max_z = np.amax(var_x), np.amax(var_y), np.amax(var_z)
            min_x, min_y, min_z = np.amin(var_x), np.amin(var_y), np.amin(var_z)
            max_i = np.amax(np.array((max_x,max_y,max_z)))
            min_i = np.amin(np.array((min_x,min_y,min_z)))
            if max_i >= max:
                    max = max_i
            if min_i <= min:
                    min = min_i
    
            #Interpolate section for all y values
            for i in range(nx):
                    for k in range(nz):
                            var_y = var[i,:,k]
                            y = np.arange(0,ny) # Arrange y values in array
                            f = interpolate.interp1d(y,var_y) # Interpolate the data
                            y_new = np.arange(0,(ny-1),step) # Arrange new y values 
                            var_y_new = f(y_new) # interpolate y values
                            var_new[i,:,k] = var_y_new # Store values in new variable
    
            #Convert coordinate section             
            vrbl = visual.vtk_var(var_new,nx,num_new,nz) # vtk variable
            vtk_path = visual.write_vtk(name,pts,vrbl,q,nx,num_new,nz) # write vtk file
            print "At t = %d, %d Steps remaining" % (q,t-(q+1)) # Progress indicator
    
    #Write the Max and min values to file
    mm_array = np.array((max,min))
    np.savetxt('max_min_' + name + '.txt',mm_array)