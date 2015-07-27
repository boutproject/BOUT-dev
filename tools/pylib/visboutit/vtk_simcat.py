#!/usr/bin/python

'''
vtk_simcat.py
Inputs:
1 - Desired end time of data to plot, if this is greater than max_time available set to max
2 - name of variable in data to plot

Outputs:
a batch of .vtk files stored in wkd/batch which is used for drawing the .vtk files
'''
#Import all the libraries
import visual
import sys
import os
import circle
import numpy as np
from scipy import interpolate
import settings as set

# Get the name of the variable to import
name = str(sys.argv[1])

#Get the working dir
work_dir = os.getcwd()

#Set the working dir to the Parent of the VisBOUTIt folder
os.chdir(work_dir +'/../')
#Update the working dir
work_dir = os.getcwd()


# Find the max time value
max_t = visual.time_max(name)

# Import coordinates from gridfile
#make this so it can find a grd file?
r = visual.nc_var(set.grid_file,'Rxy') # toroidal coordinates
z = visual.nc_var(set.grid_file,'Zxy') # toroidal coordinates

#Find Dimensions of variable data
t = 0
var_0 = visual.var3d(name,t) # collect variable
nx, ny, nz = len(var_0[:,0,0]), len(var_0[0,:,0]), len(var_0[0,0,:]) # size of grid
pts = visual.cylinder(r, z, nx, ny, nz)  # vtk grid points

#Interpolate setup
step = 0.5 #interval for new points
#Get number of new points
num_new = len(np.arange(0,(ny-1),step))
#Define empty array with increased number of y values
var_new = np.empty((nx,num_new,nz),dtype=float)

#Get end time from user input
#t = int(raw_input('Enter end time value (Enter -1 for the entire dataset): '))
t = int(sys.argv[2])

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


'''

i = 0 
#For the entire t range import the spacial values, coordinate transform, write to vtk and plot the points
for i in range(t):
	var = visual.var3d(name,i) # collect variable
	vrbl = visual.vtk_var(var,nx,ny,nz) # vtk variable
	visual.write_vtk(name,pts,vrbl,i,nx,ny,nz) # write vtk file

'''
