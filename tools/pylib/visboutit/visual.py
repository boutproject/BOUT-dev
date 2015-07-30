'''
visual.py
This file contains a library of functions used commonly by the various Scripts
'''

#==============================================================================
# Import section
#==============================================================================

from boututils import DataFile #Does this need to be imported?
from boutdata import collect

from numpy import shape #Get rid of this import?
from scipy.io import netcdf
from math import sin, cos, pi
import numpy as np
import sys
import os
from scipy import interpolate #Does this need to be imported?

#Import the  evtk library
from evtk.hl import gridToVTK
#linux Version

#==============================================================================
# # Import the vtk library mac version
# from pyevtk.hl import gridToVTK
#==============================================================================


#Import settings
import ConfigParser as cp

# Read variables from the setup file
bout_path = os.path.expanduser('~') + '/BOUT-dev' #Initial BOUT dir
if os.path.exists(bout_path):
    bout_path = bout_path
if not os.path.exists(bout_path):
    bout_path = str(raw_input('\n BOUT-dev folder not found, please enter path to BOUT-dev: '))

visboutit_path = bout_path + '/tools/pylib/visboutit/'
parser = cp.ConfigParser()
parser.read(visboutit_path + "visit.ini")
visit_dir = parser.get("file_locations","visit_dir")
visit_bin = parser.get("file_locations","visit_bin")

img_height = int(parser.get("image_settings",'img_height'))
img_width = int(parser.get("image_settings",'img_width'))


# Import the VisIt library
sys.path.insert(0,visit_dir)
import visit

#==============================================================================
# Start of the Functions
#==============================================================================

# Get returns the grid file name
    
def get(filename, name, section=None):
    
    with open(filename, "rt") as f:
        if section is not None:
            # First find the section
            found = False
            for line in f:
                # Strip spaces from left
                line = line.lstrip(' \t\n\r')
                if len(line) < 1:
                    continue  # Empty line
                    
                # if line starts with '[' then this is a section
                if line[0] == '[':
                    # Split on ']'
                    head, _ = line[1:].split(']', 1)
                    # head is now the section name
                    if head == section:
                        found = True
                        break
            if not found:
                raise ValueError("Section '%s' not found" % (section))
        
        # Now in the correct section
        
        for line in f:
            # Strip spaces from left
            line = line.lstrip(' \t\n\r')
            if len(line) < 1:
                continue  # Empty line
                
            # if line starts with '[' then this is a section
            if line[0] == '[':
                raise ValueError("Name '%s' not found in section '%s'" % (name,section))
            # Check if this line contains an '='
            if '=' in line:
                # Check if contains comment
                comment = ''
                if '#' in line:
                    line, comment = line.split('#', 1)
                # Split on '='
                key, value = line.split('=',1)
                # Strip whitespace
                key   = key.strip(' \t\n\r')
                value = value.strip(' \t\n\r')
                
                if '.' in line:
                    value, extension = value.split('.',1)
                
                # Strip out quotes if present
                if value[0] == '"' or value[0] == "'": 
                    value = value[1:]
                if value[-1] == '"' or value[-1] == "'":
                    value = value[:-1]
                
                #print("'%s' = '%s'" % (key, value))
                if key.lower() == name.lower(): # Case insensitive
                    return value


# collect a 4d variable
def var4d(name):
        var = collect(name)[:]
        return var

# NetCDF format

# collect a variable from particular .nc file
def nc_var(fname,vname):
    file = netcdf.netcdf_file(fname, 'r')
    var = file.variables[vname][:]
    file.close()
    return var

# create cylinder in .nc file
def cylinder_nc(r, nx, ny, nz):
    x_coord = []
    y_coord = []
    for i in xrange(nx):
        for k in xrange(nz):
            angle = (2./3.) * pi * (float(k)/float(nz-1))
            for j in xrange(ny):
                x_coord.append(r[i,j] * cos(angle))
                y_coord.append(r[i,j] * sin(angle))
    return x_coord, y_coord

# create torus-cylinder in .nc file
def torus_nc(r, nx, ny, nz):
    x_coord = []
    y_coord = []
    for i in xrange(nx):
        for k in xrange(nz):
            angle = (2./3.) * pi * (float(k)/float(nz-1))
            for j in xrange(ny):
                x_coord.append(r[i,j] * cos(angle))
                y_coord.append(r[i,j] * sin(angle))
    return x_coord, y_coord

# write the new .nc file
def write_file(name):
    file = netcdf.netcdf_file(name, 'w')
    return file

# write the dimension and coordinates
def coords(file, x_coord, y_coord, nt, nx, ny, nz):
        # definition of names and size of the axis
    dimnames=('level','longitude','latitude')
    dimsizes=(ny,nz,nx)
        # creation of the dims
    file.createDimension('level', ny)
    file.createDimension('longitude', nz)
    file.createDimension('latitude', nx)
    file.createDimension('time', nt)
        # creation of points in grid
    x = file.createVariable('LONGXY', 'd', dimnames)
    y = file.createVariable('LATIXY', 'd', dimnames)
    x[:] = np.reshape(x_coord,dimsizes,'F')
    y[:] = np.reshape(y_coord,dimsizes,'F')

# write the new variable
def write_var(file, var, name, nt, nx, ny, nz):
    dimnames=('time','level','longitude','latitude')
    dimsizes=(nt,ny,nz,nx)
        # necessary swap to match grid axis (if necessary)
    var = np.swapaxes(var,1,2)
    var = np.swapaxes(var,2,3)
        # write the new var
    data = file.createVariable(name, 'd', dimnames)
    data[:] = np.reshape(var,dimsizes,'F')

### VTK format



# Interpolation of 2D arrays (i.e. Rxy,Zxy) taking into account zshift
# and performing irregular number of inserts
## Inputs:
# nx,ny: number of (x,y)  points in original data
# zshift: Imported from grid file, shift in z for each (x,y) point
# z_tol: z_shift tolerance
# var: variable to interpolate
### Outputs:
# var2:  variable that has been interpolated (using linear Interpolation)
# ny2: new length of the y array

def zshift_interp2d(nx,ny,zshift,z_tol, var):
	# Input array Rxy[nx, ny]
	# Decide how many points to insert between point y and y+1
	# using zShift[nx,ny]
	# ninsert is an array of length ny, type int
	ninsert = np.zeros((ny,), dtype=np.int)
	for y in range(ny):
		if y == (ny-1):
			ninsert[y] = 0
		if y <  (ny-1):
			ninsert[y] = int( np.amax(np.abs(zshift[:,y+1] - zshift[:,y])) / z_tol )

	# Total number of points to insert
	nadd = np.sum(ninsert)

	ny2 = ny + nadd # New size of array
	var2 = np.zeros( (nx, ny2) ) # New array to put interpolated values into
	# Go through original y points
	y2 = 0  # Counter for location in new array
	
	for y in range(ny):
		# Copy y slice
		var2[:,y2] = var[:,y]
#		y2+=1
		# Interpolate and insert additional points between y and y+1
		if y < (ny-1): # Check if index is less than the length of the array
			for i in range(ninsert[y]+1):
#				var2[:,y2] = var[:,y]
				# Interpolation weights
				a = ((float(i+1)) / (ninsert[y] + 1))
				b = float((ninsert[y] - i)) / float((ninsert[y] + 1))
				y2+=1
				var2[:,y2] = (a * var[:,y+1])  + (b * var[:,y])
	return var2, ny2


# zshift_interp3d funct
def zshift_interp3d(nx ,ny ,nz ,zshift ,z_tol, var):
	# Input array var[nx,ny,nz]
	# Determine how many points to interpret between y and y+1
	# Using zshift[nx,ny]
	
	ninsert = np.zeros((ny,) , dtype=np.int)
	for y in range(ny):
		if y == (ny-1):
			ninsert[y]=0
		if y < (ny-1):
			ninsert[y] = int( np.amax(np.abs(zshift[:,y+1] - zshift[:,y]))/z_tol)

	# Total number of points to insert
	nadd = np.sum(ninsert)
	ny2 = ny + nadd # New size of array
	var2 = np.zeros( (nx, ny2, nz) ) # New array to put interpolated values into
        # Go through original y points
        y2 = 0  # Counter for location in new array

	for y in range(ny):
		# Copy y slice
		var2[:,y2,:] = var[:,y,:]
#            	   y2+=1
		# Interpolate an'd insert additional points between y and y+1
		if y < (ny-1): # Check if index is less than the length of the array
			for i in range(ninsert[y]+1):
				#var2[:,y2] = var[:,y]
				# Interpolation weights
				a = ((float(i+1)) / (ninsert[y] + 1))
				b = float((ninsert[y] - i)) / float((ninsert[y] + 1))
				y2+=1
				var2[:,y2,:] = (a * var[:,y+1,:])  + (b * var[:,y,:])
        return var2, ny2



# Creates cylinder mesh from coordinates r,z and size of grid
# Returns x,y,z coordinates as np arrays
def cylinder(r, z, nx, ny, nz):
        xcoord = np.empty((nx,ny,nz),dtype=float)
        ycoord = np.empty((nx,ny,nz),dtype=float)
        zcoord = np.empty((nx,ny,nz),dtype=float)
        # Transform Rxyz, Zxyz coordinates to cartesian
        for i in range(nx): # latitude
                for k in range(nz): # longtitude
                        phi = (2./3.)*pi*(float(k)/float(nz-1)) # (2/3)pi cylinder revolution
                        for j in range(ny): # level
                                xcoord[i,j,k] = float(r[i,j]*cos(phi))
                                ycoord[i,j,k] = float(r[i,j]*sin(phi))
                                zcoord[i,j,k] = float(z[i,j])
                                
        return xcoord,ycoord,zcoord

# Creates torus mesh from coordinates r,z and size of grid
# Returns x,y,z coordinates as np arrays
def torus(r, z, nx, ny, nz):
	xcoord = np.empty((nx,ny,nz), dtype=float)
	ycoord = np.empty((nx,ny,nz), dtype=float)
	zcoord = np.empty((nx,ny,nz), dtype=float)
	# Transform Rxyz, Zxyz coordinates to cartesian
	for i in range(nx): # latitude
    		for k in range(nz): # longtitude
        		for j in range(ny): # level
				phi = 2.*pi*(float(j)/float(ny-1)) # 2pi torus revolution
# HAVE I MESSED THIS UP A BIT WITH ORDER OF INDICIES?!!?
				xcoord[i,j,k] = float(r[i,j,k]*cos(phi))
	    			ycoord[i,j,k] = float(r[i,j,k]*sin(phi))
				zcoord[i,j,k] = float(z[i,j,k])
	return xcoord,ycoord,zcoord

#Creates ELM mesh from coordinates r,z, and size of grid
#Input Rxy,Zxy,zShift from gridfile, size of grid, and periodicity
#Output x,y,z coordinates
def elm(Rxy, Zxy, zshift, nx, ny, nz, period=1):
    dz = 2.*pi / (period*(nz-1)) # Change in z
    phi0 = np.linspace(0,2.*pi / period, nz) # Initial Phi_0 values?

	#Create empty x,y,z coordinate arraays
    xcoord = np.empty((nx,ny,nz), dtype=float)
    ycoord = np.empty((nx,ny,nz), dtype=float)
    zcoord = np.empty((nx,ny,nz), dtype=float)
	#Assign the points
#	start = 0
#	for y_i in range(ny):
#		end = start + nx*nz
#		phi = zshift[:,y_i] + phi0[:,None]
#	        r = Rxy[:,y_i] + (np.zeros([nz]))[:,None]
#	        xz_points = points[start:end]
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                #phi = zshift[i,j] + phi0[j]
                phi = (k * dz) + zshift[i,j]
                r = Rxy[i,j]
                xcoord[i,j,k] = (r*cos(phi))  # X
                ycoord[i,j,k] = (r*sin(phi))  # Y
                zcoord[i,j,k] = (Zxy[i,j])    # Z
    return xcoord,ycoord,zcoord


#Find the maximum and minimum zshift values from zshift file
def z_shift_mm(nx,ny,zshift):
	#Find max z shift
	max_z = 0
	for i in range (nx):
		for j in range(ny-1):
			z_shift = abs(zshift[i,j] - zshift[i,j+1]) #Calc zshift diff
			if z_shift > max_z:
				max_z = z_shift #Assign Max shift
	#Find min zshift
	min_z = max_z
	for i in range (nx): #For all x
		for j in range(ny-1): #For all y
			z_shift = abs(zshift[i,j] - zshift[i,j+1]) #Calc zshift diff
			if z_shift < min_z:
				min_z = z_shift #Assign shift value
	return max_z, min_z

#Find the z_tol value
#tol is a percentage of maximum z_shift is the tolerance value
#z_tol is absolute value of tolerance (above this value will get interpolated, below ignored)
def z_shift_tol(nx,ny,zshift,tol):
	#Find max z shift
        max_z = 0
        for i in range (nx):
                for j in range(ny-1):
                        z_shift = abs(zshift[i,j] - zshift[i,j+1]) #Calc zshift diff
                        if z_shift > max_z:
                                max_z = z_shift #Assign Max shift
        #Find min zshift
        min_z = max_z
        for i in range (nx): #For all x
                for j in range(ny-1): #For all y
                        z_shift = abs(zshift[i,j] - zshift[i,j+1]) #Calc zshift diff
                        if z_shift < min_z:
                                min_z = z_shift #Assign shift value
	z_tol = min_z + (max_z *tol)
	return z_tol


# Return spacial part of 4d specified variable
def var3d(name,t):
	var = collect(name)[t,:]
	return var

# Collect specified variable and return number of time slices
def time_max(name):
	var = collect(name)
	max_t = var.shape[0]
	return max_t

#Find the Dimensions of the variable data
#input: name,
def dimd(name):
    var = collect(name)
    var_0 = var[0,:]
    nx,ny,nz = len(var_0[:,0,0]), len(var_0[0,:,0]), len(var_0[0,0,:])
    return nx,ny,nz

# Find the maximum time of data, the initial value and shape (nx,ny,nz)
def dim_all(name):
    var = collect(name)
    max_t = var.shape[0]
    var_0 = var[0,:]
    nx,ny,nz = len(var_0[:,0,0]), len(var_0[0,:,0]), len(var_0[0,0,:])
    return max_t, var_0, nx, ny, nz


# create the vtk variable
def vtk_var(var, nx, ny, nz):
    vrbl = np.empty((nx,ny,nz),dtype=float)
    for i in range(nx): # latitude
        for k in range(nz): # longtitude
            for j in range(ny): # level
                vrbl[i,j,k] = var[i,j,k] #i,j,k or i,k,j?
    return vrbl

# Writes data to vtk file for every time slice
# Returns the vtk file path
def write_vtk(name,pts,vrbl,t,nx,ny,nz):
	i,j,k = pts
	vrbl = np.array(vrbl)
	vtk_file_path = gridToVTK("./batch/" + name + "_batch_%d" % t, i, j, k, pointData = {str(name) : vrbl})
	return vtk_file_path



# Write the real and imaginary data to vtk file
# Inputs:
# name: name of variable
# pts: grid points
# vrbl_r: real eigenvalues
# vrbl_i: imaginary eigenvalues
# eig_num: eigen number,
def write_vtk_2(name,pts,vrbl_r,vrbl_i,eig_num):
	i,j,k = pts
#	names = (name + '_r', name + '_i')
#	vrbl_m = vrbl_r,vrbl_i
	vtk_file_path = gridToVTK( name + "_eigen_%d" % eig_num,i,j,k, pointData = {str(name + '_r'): vrbl_r , str(name + '_i') : vrbl_i})
	return vtk_file_path



#Draw vtk file and let user orientate view and then save session file
# returns the VisIt session name and location
def view_vtk(work_dir,name,max,min):
    vtk_path = work_dir + "/batch/" + name + "_batch_*.vts database" #Set vtkfile path
    visit.OpenDatabase(vtk_path) # Open database
    visit.AddPlot("Pseudocolor",name) #Draw a Pseudocolor Plot of the variable
    # If user would like fixed max and min then assign max and min
    #Set the max and min values for the data
    PseudocolorAtts = visit.PseudocolorAttributes()
    if max != False:
        PseudocolorAtts.max = max
        PseudocolorAtts.maxFlag = 1
        visit.SetPlotOptions(PseudocolorAtts)
    if min != False:
        PseudocolorAtts.min = min
        PseudocolorAtts.minFlag = 1
        visit.SetPlotOptions(PseudocolorAtts)
        visit.DrawPlots() #Draw the plots
    #Save the Visit Session
    session_name = raw_input('Enter a session file name:')
    session_path = work_dir+ "/" + session_name + ".session"
    visit.SaveSession(session_path)
    #Close visit session
    visit.DeleteAllPlots()
    visit.CloseDatabase(vtk_path)
    return session_path,session_name

#Export an image sequence of the plot across the entire time range
def draw_vtk(session_path,img_dir,name,t,session_name,max_imp,min_imp):
    #Make dir for storing image sequence
    outputdir = img_dir + '/' + session_name
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if max_imp == False:
        max_imp = 0
    if min_imp == False:
        min_imp = 0

    #Launch visit
    sys.path.insert(0,visit_dir)
    import visit
    #Load session and initialise at time 0
    visit.RestoreSession(session_path,0)
    visit.SetTimeSliderState(0)
    #Export an Image sequence of the variable for every time base
    i = 0
    for i in range(t):
        visit.SetTimeSliderState(i) #Change timer slider
        #Make this more general to different plots!?
        PseudocolorAtts = visit.PseudocolorAttributes()
        # If user would like fixed max and mind then assign max and min values
        if max_imp != 0:
            PseudocolorAtts.max = max_imp
            PseudocolorAtts.maxFlag = 1
        if max_imp == 0:
            PseudocolorAtts.maxFlag = 0
        if min_imp != 0:
            PseudocolorAtts.min = min_imp
            PseudocolorAtts.minFlag = 1
        if min_imp == 0:
            PseudocolorAtts.minFlag = 0
        visit.SetPlotOptions(PseudocolorAtts)
        visit.DrawPlots() # Draw plot
        # Save a png of the plot
        s = visit.SaveWindowAttributes()
        s.outputToCurrentDirectory = 0
        s.outputDirectory = outputdir
        s.family = 0
        s.fileName = '%s_%s_image_%04d' % (name,session_name,i)
        s.format = s.PNG
        s.width = img_width
        s.height = img_height
        visit.SetSaveWindowAttributes(s)
        visit.SaveWindow()
    #Close visit session
    visit.DeleteAllPlots()
    visit.Close()