"""
visual.py
This file contains a library of functions used by various scripts in the VisBOUTIt package
"""

#==============================================================================
# Import section
#==============================================================================

from boututils.datafile import DataFile
from boutdata import collect
from scipy.io import netcdf
from numpy import sin, cos, pi
import numpy as np
import sys
import os
from scipy import interpolate


#Import the  evtk library
from evtk.hl import gridToVTK

#Import settings
import configparser as cp

#==============================================================================
# Read setup file section
#==============================================================================

# Find the BOUT-dev/tools/pylib folder from the users bash profile.
visboutit_path = None
pypaths = sys.path
for line in pypaths:
    if 'BOUT' in line:
        if 'tools' in line:
            if 'pylib' in line:
                visboutit_path = line + '/visboutit/'
    
if visboutit_path == None:
    bout_path = str(raw_input('\n BOUT-dev folder not found, please enter path to BOUT-dev: '))
    visboutit_path = bout_path + '/tools/pylib/visboutit/'

# Read variables from the setup file
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

def var_list():
    d = DataFile('BOUT.dmp.0.nc')
    var_list = d.list()
    return var_list


def zShf_p_check(zShf_int_p):
    while zShf_int_p > 1:
        print('Error: zShift_int_percent has to be positive and between 0 and 1')
        zShf_int_p = float(raw_input('Enter New zShift_int_percent value: '))
    return zShf_int_p

# collect a 4d variable
def var4d(name):
        var = collect(name)[:]
        return var

#==============================================================================
#  NetCDF format
#==============================================================================

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

#==============================================================================
#  VTK format
#==============================================================================

def intrp_grd(nx, ny, var, ny_work, step):
    if len(var.shape) == 2:
        var_new = np.empty((nx, ny_work),dtype = float)
    elif len(var.shape) == 3:
        nz = var.shape[2]
        var_new = np.empty((nx, ny_work, nz),dtype = float)
    else:
        raise ValueError("intrp_grd expects 2d or 3d variable")

    for i in range(nx):
        var_yslice = var[i,:]
        y = np.arange(0,ny)
        f = interpolate.interp1d(y, var_yslice, axis=0) # axis zero because we restricted to x-index i already
        y_new = np.arange(0, (ny-1), step)
        var_y_new = f(y_new)

        var_new[i,:] = var_y_new
    return var_new


#Interpolation of 2D arrays (i.e. Rxy,Zxy) taking into account zshift
#and performing irregular number of inserts
#Inputs:
#nx,ny: number of (x,y)  points in original data
#zshift: Imported from grid file, shift in z for each (x,y) point
#z_tol: z_shift tolerance
#var: variable to interpolate
#Outputs:
#var2:  variable that has been interpolated (using linear Interpolation)
#ny2: new length of the y array


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
        # Interpolate and insert additional points between y and y+1
        if y < (ny-1): # Check if index is less than the length of the array
            for i in range(ninsert[y]+1):
                # Interpolation weights
                a = ((float(i+1)) / (ninsert[y] + 1))
                b = float((ninsert[y] - i)) / float((ninsert[y] + 1))
                y2+=1
                var2[:,y2] = (a * var[:,y+1])  + (b * var[:,y])
    return var2, ny2

# zshift_interp3d funct
def zshift_interp3d(nx ,ny ,nz ,zshift ,z_tol, var):
    """
    Input array var[nx,ny,nz]
    Determine how many points to interpret between y and y+1
    Using zshift[nx,ny]
    """
    
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
        # Interpolate an'd insert additional points between y and y+1
        if y < (ny-1): # Check if index is less than the length of the array
            for i in range(ninsert[y]+1):
                # Interpolation weights
                a = ((float(i+1)) / (ninsert[y] + 1))
                b = float((ninsert[y] - i)) / float((ninsert[y] + 1))
                y2+=1
                var2[:,y2,:] = (a * var[:,y+1,:])  + (b * var[:,y,:])
    return var2, ny2

# Creates slab mesh from coordinates r,z and size of grid
# Returns x,y,z coordinates as np arrays
def slab(x, y, zShift, nx, ny, nz):
    xcoord = np.empty((nx,ny,nz),dtype=float)
    ycoord = np.empty((nx,ny,nz),dtype=float)
    zcoord = np.empty((nx,ny,nz),dtype=float)
    z = 
    # Transform Rxyz, Zxyz coordinates to cartesian
    xcoord[...] = x[:, :, np.newaxis]
    ycoord[...] = y[:, :, np.newaxis]
    zcoord[...] = z[np.newaxis, np.newaxis, :] + zShift[:, :, np.newaxis]
                                
    return xcoord,ycoord,zcoord

# Creates cylinder mesh from coordinates r,z and size of grid
# Returns x,y,z coordinates as np arrays
def cylinder(r, z, nx, ny, nz, pi_fr = (2./3.)):
    # pi_fr check if greater than 2
    while pi_fr > 2:
        pi_fr = pi_fr - 2
    
    xcoord = np.empty((nx,ny,nz),dtype=float)
    ycoord = np.empty((nx,ny,nz),dtype=float)
    zcoord = np.empty((nx,ny,nz),dtype=float)
    pi_fr = float(pi_fr)
    phi = np.linspace(0., pi_fr*pi/nz, nz, endpoint=False)[np.newaxis, np.newaxis, :]
    # Transform Rxyz, Zxyz coordinates to cartesian
    xcoord[...] = r[:, :, np.newaxis] * cos(phi)
    ycoord[...] = r[:, :, np.newaxis] * sin(phi)
    zcoord[...] = z[:, :, np.newaxis]
                                
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
                xcoord[i,j,k] = float(r[i,j,k]*cos(phi))
                ycoord[i,j,k] = float(r[i,j,k]*sin(phi))
                zcoord[i,j,k] = float(z[i,j,k])
    return xcoord,ycoord,zcoord

# Creates ELM mesh from coordinates r,z, and size of grid
# Input Rxy,Zxy,zShift from gridfile, size of grid, and periodicity
# Output x,y,z coordinates
def elm(Rxy, Zxy, zshift, nx, ny, nz, period=1):
    dz = 2.*pi / (period*(nz-1)) # Change in z
    phi0 = np.linspace(0,2.*pi / period, nz) # Initial Phi_0 value
    # Create empty x,y,z coordinate arraays
    xcoord = np.empty((nx,ny,nz), dtype=float)
    ycoord = np.empty((nx,ny,nz), dtype=float)
    zcoord = np.empty((nx,ny,nz), dtype=float)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                phi = (k * dz) + zshift[i,j]
                r = Rxy[i,j]
                xcoord[i,j,k] = (r*cos(phi))  # X
                ycoord[i,j,k] = (r*sin(phi))  # Y
                zcoord[i,j,k] = (Zxy[i,j])    # Z
    return xcoord,ycoord,zcoord

# Find the maximum and minimum zshift values from zshift file
def z_shift_mm(nx,ny,zshift):
    # Find max z shift
    max_z = 0
    for i in range (nx):
        for j in range(ny-1):
            z_shift = abs(zshift[i,j] - zshift[i,j+1]) # Calc zshift diff
            if z_shift > max_z:
                max_z = z_shift # Assign Max shift
    # Find min zshift
    min_z = max_z
    for i in range (nx): # For all x
        for j in range(ny-1): # For all y
            z_shift = abs(zshift[i,j] - zshift[i,j+1]) # Calc zshift diff
            if z_shift < min_z:
                min_z = z_shift # Assign shift value
    return max_z, min_z

# Find the z_tol value
# tol is a percentage of maximum z_shift is the tolerance value
# z_tol is absolute value of tolerance (above this value will get interpolated, below ignored)
def z_shift_tol(nx,ny,zshift,tol):
    # Find max z shift
    max_z = 0
    for i in range (nx):
        for j in range(ny-1):
            z_shift = abs(zshift[i,j] - zshift[i,j+1]) # Calc zshift diff
            if z_shift > max_z:
                max_z = z_shift # Assign Max shift
    # Find min zshift
    min_z = max_z
    for i in range (nx): # For all x
        for j in range(ny-1): # For all y
            z_shift = abs(zshift[i,j] - zshift[i,j+1]) # Calc zshift diff
            if z_shift < min_z:
                min_z = z_shift # Assign shift value
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

# Find the Dimensions of the variable data
# input: name
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
    nx,ny,nz = var_0.shape
    return max_t, var_0, nx, ny, nz

# create the vtk variable
def vtk_var(var, nx, ny, nz):
    vrbl = np.empty((nx,ny,nz),dtype=float)
    for i in range(nx): # latitude
        for k in range(nz): # longtitude
            for j in range(ny): # level
                vrbl[i,j,k] = var[i,j,k]
    return vrbl

# Writes data to vtk file for every time slice
# Returns the vtk file path
def write_vtk(name,pts,vrbl,t):
    i,j,k = pts
    vrbl = np.array(vrbl)
    vtk_file_path = gridToVTK("./batch/" + name + "_batch_%d" % t, i, j, k, pointData = {str(name) : vrbl})
    return vtk_file_path

# Write the real and imaginary data to vtk file
# Inputs:
# name: name of variables
# pts: grid points
# vrbl_r: real eigenvalues
# vrbl_i: imaginary eigenvalues
# eig_num: eigen number,
def write_vtk_2(name,pts,vrbl_r,vrbl_i,eig_num):
    i,j,k = pts
    vtk_file_path = gridToVTK("./batch/" + name + "_eigen_%d" % eig_num,i,j,k, pointData = {str(name + '_r'): vrbl_r , str(name + '_i') : vrbl_i})
    return vtk_file_path

def write_vtk_vector(name,pts,vectorx,vectory,vectorz,t):
    i,j,k = pts
    vtk_file_path = gridToVTK("./batch/" + name + "_vector_%d" % t , i , j , k , pointData = {str(name): (vectorx,vectory,vectorz)})
    return vtk_file_path

def write_vtk_vector2(name,pts,vector,vectorx,vectory,vectorz,t):
    i,j,k = pts
    vtk_file_path = gridToVTK("./batch/" + name + "_vector_%d" % t , i , j , k , pointData = {str(name):vector, str(name + '_vector'): (vectorx,vectory,vectorz)})
    return vtk_file_path

def write_vtk_vector3(name, pts, vectora, vectorb, vectorc):
    i,j,k = pts
    vectora_i, vectora_j, vectora_k = vectora
    vectorb_i, vectorb_j, vectorb_k = vectorb
    vectorc_i, vectorc_j, vectorc_k = vectorc
    vtk_file_path = gridToVTK("./batch/" + name + "_vector" , i , j , k , pointData = {str(name+ '_a'): (vectora_i, vectora_j, vectora_k), str(name + '_b'): (vectorb_i, vectorb_j, vectorb_k), str(name + '_c'): (vectorc_i, vectorc_j, vectorc_k)})
    return vtk_file_path

# Draw vtk file and let user orientate view and then save session file
# returns the VisIt session name and location
def view_vtk(work_dir,name,max,min):
    vtk_path = work_dir + "/batch/" + name + "_batch_*.vts database" # Set vtkfile path
    visit.OpenDatabase(vtk_path) # Open database
    visit.AddPlot("Pseudocolor",name) #Draw a Pseudocolor Plot of the variable
    # If user would like fixed max and min then assign max and min
    # Set the max and min values for the data
    PseudocolorAtts = visit.PseudocolorAttributes()
    if max != False:
        PseudocolorAtts.max = max
        PseudocolorAtts.maxFlag = 1
        visit.SetPlotOptions(PseudocolorAtts)
    if min != False:
        PseudocolorAtts.min = min
        PseudocolorAtts.minFlag = 1
        visit.SetPlotOptions(PseudocolorAtts)
        visit.DrawPlots() # Draw the plots
    # Save the Visit Session
    session_name = raw_input('Enter a session file name:')
    session_path = work_dir+ "/" + session_name + ".session"
    visit.SaveSession(session_path)
    # Close visit session
    visit.DeleteAllPlots()
    visit.CloseDatabase(vtk_path)
    return session_path,session_name

# Export an image sequence of the plot across the entire time range
def draw_vtk(session_path,img_dir,name,t,session_name,max_imp,min_imp,skip):
    # Make dir for storing image sequence
    outputdir = img_dir + '/' + session_name
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if max_imp == False:
        max_imp = 0
    if min_imp == False:
        min_imp = 0

    # Create index of number of time slices (used for keeping time in the filenames and navigating in VisIt)
    indicies = np.arange(0,len(t))
        
    # Launch visit
    sys.path.insert(0,visit_dir)
    import visit
    # Load session and initialise at time 0
    visit.RestoreSession(session_path,0)
    visit.SetTimeSliderState(0)
    # Export an Image sequence of the variable for every time base
    for i in indicies:
        time = i * skip
        visit.SetTimeSliderState(i) # Change timer slider
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
        s.fileName = '%s_%s_image_%04d' % (name,session_name,time)
        s.format = s.PNG
        s.width = img_width
        s.height = img_height
        visit.SetSaveWindowAttributes(s)
        visit.SaveWindow()
    # Close visit session
    visit.DeleteAllPlots()
    visit.Close()
    
# Draw vtk file and let user orientate view and then save session file
# returns the VisIt session name and location
def view_eigen(work_dir,name):
    vtk_path = work_dir + "/batch/" + name + "_eigen_*.vts database" # Set vtkfile path
    visit.OpenDatabase(vtk_path) # Open database
    # Create two windows for real and imaginary values
    visit.SetWindowLayout(2)
    visit.SetActiveWindow(1)
    visit.AddPlot("Pseudocolor",name + '_r', 1 , 1) #Draw a Pseudocolor Plot of the real eigen variable
    visit.DrawPlots() #Draw real plot
    visit.SetActiveWindow(2)
    visit.AddPlot("Pseudocolor",name + '_i', 1 , 1) #Draw a Pseudocolor Plot of the imaginary eigen variable
    visit.DrawPlots() #Draw imaginary plot       
        
    # Save the Visit Session
    session_name = raw_input('Enter a session file name:')
    
    #Detect if 2D or 3D 
    visit.SetActiveWindow(1)
    dim_1 = visit.GetWindowInformation().viewDimension
    visit.SetActiveWindow(2)
    dim_2 = visit.GetWindowInformation().viewDimension
    
    # Return error if dimensions of the plots are not the same
    if dim_1 != dim_2:
        print('Error plots do not have the same dimensions')
        sys.exit()
    
    # Set dimension variable
    if dim_1 == dim_2:
        dim = dim_1
        
    # Set both windows to have the same orientation as window 1, 2D
    if dim == 2:
        visit.SetActiveWindow(1)
        view = visit.View2DAttributes(1)
        visit.SetActiveWindow(2)
        visit.SetVeiew2D(view)
    

    # Set both windows to have the same orientation as window 1, 3D
    if dim == 3: 
        visit.SetActiveWindow(1)
        view = visit.View3DAttributes(1)
        visit.SetActiveWindow(2)
        visit.SetView3D(view)
    
    # Save the session
    session_path = work_dir + "/" + session_name + ".session"
    visit.SaveSession(session_path)
    
    # Delete plots
    for i in (1,2):
        visit.SetActiveWindow(i)
        visit.DeleteAllPlots()
    
    visit.CloseDatabase(vtk_path)
    return session_path,session_name
    
    
# Function that renders an image sequence of a session with two windows (designed for eigen values)
#==============================================================================
# Create an image sequence of eigen data
#==============================================================================
def image_eigen(session_path,img_dir,name,t,session_name):
    # Make dir for storing image sequence
    outputdir = img_dir + '/' + session_name
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)


    # Create index of number of time slices (used for keeping time in the filenames and navigating in VisIt)
    index = np.arange(0,len(t))
        
    # Launch visit
    sys.path.insert(0,visit_dir)
    import visit
    # Load session and initialise at time 0
    visit.RestoreSession(session_path,0)
    visit.SetTimeSliderState(0)
    
    # Export an Image sequence of the real variable for every time base
    for i in index:
        visit.SetActiveWindow(1)
        visit.SetTimeSliderState(i) # Change timer slider
        visit.DrawPlots() # Draw plot
        # Save a png of the plot
        s = visit.SaveWindowAttributes()
        s.outputToCurrentDirectory = 0
        s.outputDirectory = outputdir
        s.family = 0
        s.fileName = '%s_%s_image_%04d' % (name + '_r',session_name,t[i])
        s.format = s.PNG
        s.width = img_width
        s.height = img_height
        visit.SetSaveWindowAttributes(s)
        visit.SaveWindow()
        
    # Export an image sequence of the real variable for every time base
    for i in index:
        visit.SetActiveWindow(2)
        visit.SetTimeSliderState(i)
        
        visit.DrawPlots()
        s = visit.SaveWindowAttributes()
        s.outputToCurrentDirectory = 0
        s.outputDirectory = outputdir
        s.family = 0
        s.fileName = '%s_%s_image_%04d' % (name + '_i',session_name,t[i])
        s.format = s.PNG
        s.width = img_width
        s.height = img_height
        visit.SetSaveWindowAttributes(s)
        visit.SaveWindow()
        
    # Close visit session
    for i in (1,2):
        visit.SetActiveWindow(i)
        visit.DeleteAllPlots()
        
    visit.Close()
