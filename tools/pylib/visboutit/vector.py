# -*- coding: utf-8 -*-
"""
vectory.py Draws vectors that are specified.
"""
#Import the relevant libraries
from visboutit import visual, circle
import os
import numpy as np
from math import sin,cos,pi,fabs
from scipy import interpolate

#==============================================================================
# Unit vector functions
#==============================================================================

def torus_y_unit(ny):
    """
    torus_y_unit returns a list of unit vectors along the field lines in a toroidal geometry.
    (clockwise model)
    
    """
    y_unit = []
    y_unit_x = np.empty((ny))
    y_unit_y = np.empty((ny))
    # Find initial y_unit
    delta_phi = 2*pi*1/(ny-1)
    y_unit_initial =  np.dot((1./(delta_phi)) , [-cos(0) + cos(delta_phi) , sin(0) - sin(delta_phi)])

#Rotate initial points to find the other remaining values    
    for j in range(ny):
        #Calculate the angle transversed
        phi = 2*pi*j/(ny-1)
        #Update the rotation matrix
        r_mtrx = [[cos(phi) , -sin(phi)] , [sin(phi) , cos(phi) ]]
        y_unit.append(np.dot(r_mtrx, y_unit_initial))
    
    return y_unit


def elm_field_vector(pts,nx,ny,nz):
    """
    elm_field_vector function returns the unit vectors along the field line of an ELM mesh
    
    """
    xcoord, ycoord, zcoord = pts
    # Create the unit vector variable
    unit_vector_x = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_y = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_z = np.zeros((nx,ny,nz),dtype = float)
    
    # For every x,y,z point calculate a unit vector along the field line
    # by calculating the vector between two points from mesh (this already accounts for zShift)
    print 'Creating vector map:'
    for i in range(nx-1):
        percent = (float(i)/(float(nx - 2)))*100
        print '%d %% Complete' % percent, "   \r", # Progress percent
        for j in range(ny-1):
            for k in range(nz-1):
                r1 = np.array((xcoord[ i , j , k ] , ycoord[ i , j , k ] , zcoord[ i , j , k ])) # Vector of initial point
                r2 = np.array((xcoord[ i , j+1 , k ] , ycoord[ i , j+1 , k ]  , zcoord[ i , j+1 , k ])) # Vector of secondary point
                r = r2-r1
                unit = r/ np.linalg.norm(r) # Normalise vector
                # Assign the vectors
                unit_vector_x[i,j,k] = unit[0]
                unit_vector_y[i,j,k] = unit[1]
                unit_vector_z[i,j,k] = unit[2]
    print "" # For closing updated print statement
    return unit_vector_x, unit_vector_y, unit_vector_z


def y_unit_vector(pts , nx , ny , nz ):
    """
    y_unit_vector function returns the y unit vectors for a given mesh
    
    Inputs:
        pts: mesh of for the coordinate system
        nx, ny ,nz: size of the mesh
        
    Output:
        unit vector in y direction in components
    
    """
    xcoord, ycoord, zcoord = pts
    # Create the unit vector variable
    unit_vector_x = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_y = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_z = np.zeros((nx,ny,nz),dtype = float)
    
    # For every x,y,z point calculate a unit vector along the field line
    # by calculating the vector between two points from mesh (this already accounts for zShift)
    print 'Creating vector map:'
    for i in range(nx-1):        
        for j in range(ny-1):
            percent = (float(i)/(float(nx - 2)))*100
            print '%d %% Complete' % percent, "   \r", # Progress percent
            for k in range(nz-1):
                r1 = np.array((xcoord[ i , j , k ] , ycoord[ i , j , k ] , zcoord[ i , j , k ])) # Vector of initial point
                r2 = np.array((xcoord[ i , j+1 , k ] , ycoord[ i , j+1 , k ]  , zcoord[ i , j+1 , k ])) # Vector of secondary point
                r = r2-r1
                unit = r/ np.linalg.norm(r) # Normalise vector
                # Assign the vectors
                unit_vector_x[i,j,k] = unit[0]
                unit_vector_y[i,j,k] = unit[1]
                unit_vector_z[i,j,k] = unit[2]
    print "" # For closing updated print statement
    return unit_vector_x, unit_vector_y, unit_vector_z


def dPhi_da(cpos , hpos , lpos, cphi, hphi, lphi):
    #dx
    dx_a = hpos[0] - cpos[0]
    dx_b = cpos[0] - lpos[0]
    dx_av = np.mean((np.array((dx_a,dx_b))))
    
#    dphi
    dphi_a = hphi - cphi
    dphi_b = cphi - lphi
    dphi_av = np.mean((np.array((dphi_a,dphi_b))))
    
    dphi_dx = dphi_av / dx_av
    
    #Get unit in x from average of r_2 -> r_0 -> r_1? or just use r_0 - > r_1?????
    rx_a = hpos - cpos
    rx_b = cpos - lpos
    rx = rx_a + rx_b
    rx_unit = rx / np.linalg.norm(rx)
    
    dphi_dx_v = dphi_dx * rx_unit
    
    return dphi_dx_v



## Function for calculating grad phi 
#def grad_phi_old(phi,pts,nx,ny,nz):
#    xcoord, ycoord, zcoord = pts
#    # Start at initial point and get the points arround it,
#    dphi_dx = []
#    dphi_dy = []
#    dphi_dz = []
#    for i in np.arange(1,nx-1):
#        print i
#        for j in np.arange(1,ny-1):
#            for k in np.arange(1,nz-1):
#                phi_0 = phi[i , j , k ]
#                phi_1 = phi[i+1 , j , k ]
#                phi_2 = phi[i-1 , j , k ]
#                phi_3 = phi[i , j , k+1 ]
#                phi_4 = phi[i , j , k-1 ]
#                phi_5 = phi[i , j+1 , k ]
#                phi_6 = phi[i , j-1 , k ]
#                r_0 = np.array((xcoord[i , j , k ] , ycoord[ i , j , k ] , zcoord[ i , j , k ])) # Initial point
#                r_1 = np.array((xcoord[i+1 , j , k ] , ycoord[ i+1 , j , k ] , zcoord[ i+1 , j , k ])) # Point one at x+1
#                r_2 = np.array((xcoord[i-1 , j , k ] , ycoord[ i-1 , j , k ] , zcoord[ i-1 , j , k ])) # Point one at x-1
#                r_3 = np.array((xcoord[i , j , k+1 ] , ycoord[ i , j , k+1 ] , zcoord[ i , j , k+1 ])) # Point one at z+1
#                r_4 = np.array((xcoord[i , j , k-1 ] , ycoord[ i , j , k-1 ] , zcoord[ i , j , k-1 ])) # Point one at z-1
#                r_5 = np.array((xcoord[i , j+1 , k ] , ycoord[ i , j+1 , k ] , zcoord[ i , j+1 , k ])) # Point one at y+1
#                r_6 = np.array((xcoord[i , j-1 , k ] , ycoord[ i , j-1 , k ] , zcoord[ i , j-1 , k ])) # Point one at y-1
#                
#
#                dphi_dx.append(dPhi_da(r_0 , r_1 , r_2 , phi_0 , phi_1 , phi_2))
#                dphi_dz.append(dPhi_da(r_0 , r_3 , r_4 , phi_0 , phi_3 , phi_4))
#                dphi_dy.append(dPhi_da(r_0 , r_5 , r_6 , phi_0 , phi_5 , phi_6))
#    return dphi_dx, dphi_dy, dphi_dz
                
                

                
# Function for calculating grad phi 
def grad_phi(phi,pts,nx,ny,nz):
    print "Grad_phi"
    # Create variables for vectors
    dphi = np.zeros((nx,ny,nz,3))
    phi_0_solve = np.zeros((nx,ny,nz))
    dphi_dx = np.zeros((nx,ny,nz))
    dphi_dy = np.zeros((nx,ny,nz))
    dphi_dz = np.zeros((nx,ny,nz))

    # Extract coordinates of mesh
    xcoord, ycoord, zcoord = pts
    
    # Start at initial point and get the points arround it,
    for i in range(nx-1):
        percent = (float(i)/(float(nx - 2)))*100
        for j in range(ny-1):
            print '%d %% Complete i=%d' % (percent,i), "   \r", # Progress percent
            for k in range(nz-1):     
                
                # Import the variables taking into accound normalisation factor and vectors needed 
                phi_0 = phi[i   , j   , k ]
                phi_1 = phi[i+1 , j   , k ]
                phi_2 = phi[i   , j+1 , k ]
                phi_3 = phi[i   , j   , k+1 ]
                r_0 = np.array((xcoord[i , j , k ]  ,  ycoord[ i , j , k ]   , zcoord[ i , j , k ])) # Initial point
                r_1 = np.array((xcoord[i+1 , j , k ] , ycoord[ i+1 , j , k ] , zcoord[ i+1 , j , k ])) # Point one at x+1
                r_2 = np.array((xcoord[i , j+1 , k ] , ycoord[ i , j+1 , k ] , zcoord[ i , j+1 , k ])) # Point one at y+1
                r_3 = np.array((xcoord[i , j , k+1 ] , ycoord[ i , j , k+1 ] , zcoord[ i , j , k+1 ])) # Point one at z+1
                
                # Solve grad_phi in components (excluding direction) using matrix
                mtrx_cnsts = np.array([[ 1 , 0 , 0 , 0 ] ,
                                       [ 1 , ( r_1[0] - r_0[0] ) , ( r_1[1] - r_0[1] ) , ( r_1[2] - r_0[2] ) ], 
                                       [ 1 , ( r_2[0] - r_0[0] ) , ( r_2[1] - r_0[1] ) , ( r_2[2] - r_0[2] ) ], 
                                       [ 1 , ( r_3[0] - r_0[0] ) , ( r_3[1] - r_0[1] ) , ( r_3[2] - r_0[2] ) ]])
                mtrx_phi = np.array([[phi_0],[phi_1],[phi_2],[phi_3]])
                mtrx_soln = np.linalg.solve(mtrx_cnsts,mtrx_phi)
                
                # Assign solutions
#                phi_0_solve[i,j,k] = mtrx_soln[0]
                dphi_dx[i,j,k] = mtrx_soln[1]
                dphi_dy[i,j,k] = mtrx_soln[2]
                dphi_dz[i,j,k] = mtrx_soln[3]
                
                #Calculate the directions of grad phi components
                i_unit_work = (r_1 - r_0)/(np.linalg.norm(r_1 - r_0))
                j_unit_work = (r_2 - r_0)/(np.linalg.norm(r_2 - r_0))
                k_unit_work = (r_3 - r_0)/(np.linalg.norm(r_3 - r_0))
                
                # calculate the dot product between maginitde and direction of the components of grad within a cartesian coordinate system
                dphi_dx_v = np.dot(dphi_dx[i,j,k] , i_unit_work)
                dphi_dy_v = np.dot(dphi_dy[i,j,k] , j_unit_work)
                dphi_dz_v = np.dot(dphi_dz[i,j,k] , k_unit_work)
                dphi_v_work = dphi_dx_v + dphi_dy_v + dphi_dz_v
                
                # Assign dphi components into dphi                
                for z in range(3):
                    dphi[i,j,k,z] = dphi_v_work[z]


    print ""
    return dphi

                

# the elm function imports

def phi_elm(name , time ,zShf_int_p = 0.25, path = None, skip = 1):
    """
       
    elm function: 
    Imports a variable that is parallel along a field line (for example; jpar) creates a mesh for
    the elm geoemtry, calculates the unit vectors along the field lines from the mesh (as this accounts for zShift).
    The function will output the vector and scalar form of this variable.
    
    
    
    Inputs
        name: variable name imported from BOUT.dmp files
        time: The end time slice that will be converted, if -1 entire dataset is imported
        zShift Interpolation Percentage: This is used to set the tolerance values for the linear interpolation.
        path: File path to data, if blank then the current working directory is used to look for data
    
    Output
        Timestamped vtk files of the specified variable, with scalar and vector variables.
        A text file with the maximum and minimum files.
                
    """
    zShf_int_p = visual.zShf_p_check(zShf_int_p)
    
    #Get the working dir
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    grid_file = visual.get("BOUT.inp","grid") # Import grid file name 
    
    # Get the dimensions and initial value of the variable data.    
    max_t, var_0, nx, ny, nz = visual.dim_all(name)
      
    #ELM
    #Import r,z and zshift from grid file
    r = visual.nc_var(grid_file,'Rxy') # ELM Coordinates
    z = visual.nc_var(grid_file,'Zxy') # ELM Coordinates
    zshift = visual.nc_var(grid_file,'zShift') #ELM Coordinates
    
    max_z,min_z = visual.z_shift_mm(nx,ny,zshift) #Find the max and minimum zshift values
    z_tol = min_z + ((max_z - min_z)* zShf_int_p) # Set tolerance value
    
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
    
#    #Initial max and min values
#    var_x, var_y, var_z = var_0[:,0,0], var_0[0,:,0], var_0[0,0,:]
#    max_x, max_y, max_z = np.amax(var_x), np.amax(var_y), np.amax(var_z)
#    min_x, min_y, min_z = np.amin(var_x), np.amin(var_y),np.amin(var_z)
#    max = np.amax((max_x, max_y, max_z))
#    min = np.amin((min_x, min_y, min_z))
    
    vectorx = np.zeros((nx,ny2,nz),dtype = float) # Create the empty vector variables
    vectory = np.zeros((nx,ny2,nz),dtype = float)
    vectorz = np.zeros((nx,ny2,nz),dtype = float)
     #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    q = 0 
    while q <= t-1:
        var = visual.var3d(name,q) # collect variable
    
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
    
    
        var2,ny3 = visual.zshift_interp3d(nx,ny,nz,zshift,z_tol,var) # Interpolate Variable
        
        vectorx , vectory , vectorz = grad_phi(var2,pts2,nx,ny2,nz)
        
        
        #Convert coordinate section
        vrbl = visual.vtk_var(var2,nx,ny2,nz) #vtk variable
        vtk_path = visual.write_vtk_vector2(name,pts2,vrbl,vectorx,vectory,vectorz,q)
        print "At t = %d, %d Steps remaining" % (q,((t- q)/skip))# Progress indicator
        q += skip
    print ""
#    #Write the Max and min values to file
#    mm_array = np.array((max,min))
#    np.savetxt('max_min_' + name + '.txt',mm_array) 
    return
                
                
                
                

def vel_torus(name, time, step = 0.5, path = None, R = None, r_min = None , dr = None , Bt = None, q = None , isttok = False, skip = 1):
    """
    Torus function writes the vector for a variable that is parallel along a field line (i.e Vpar)
    The variable is imported, interpolated, a mesh is created for the Torous geometry
    unit vectors along the field direction are created using the torus_y_unit function then using imported variable
    the vectors are created and plotted
    Output: time stamped vtk files that contain the scalar variable and a vector variable of the imported 
    
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
        Timestamped vtk files of the specified variable, scalar and vector forms.
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
    vi_new = np.empty((nx,ny_work,nz),dtype=float)
    
    if isttok == True: # Use ISTTOK specifications?
        R,r_min,dr,Bt,q = 0.46, 0.085, 0.02, 0.5, 5
    # Import coordinates from function in circle library
    r,z = circle.coordinates(nx,ny_work,nz,R,r_min,dr,Bt,q) # toroidal coordinates
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
    vectorx = np.zeros((nx,ny_work,nz),dtype= float)
    vectory = np.zeros((nx,ny_work,nz),dtype= float)
    vectorz = np.zeros((nx,ny_work,nz),dtype= float)
    
    #Find the y_unit vectors
    y_unit = torus_y_unit(ny_work)
    
    # Convert y_unit into a 3d vector
    y_unit_3d = []
    for i in y_unit:
        y_unit_3d.append(np.array((i[0],i[1],0)))


    #import normalisation factors
    t_norm = visual.collect('Tnorm')
    b_norm = visual.collect('Bnorm')

    q = 0 
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    while q <= t-1:
        var = visual.var3d(name,q)*t_norm # collect variable phi taking into account normalisation factor
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
        
        gradphi = grad_phi(var_new,pts,nx,ny_work,nz)
        
        #Import parallel velocity
        vi = visual.var3d('vi',q)
        
        #Interpolate values
        if ny != ny_work: #Interpolate section for all y values
            for i in range(nx):
                for k in range(nz):
                    vi_y = vi[i,:,k]
                    y = np.arange(0,ny) # Arrange y values in array
                    f = interpolate.interp1d(y,vi_y) # Interpolate the data
                    y_new = np.arange(0,(ny-1),step) # Arrange new y values 
                    vi_y_new = f(y_new) # interpolate y values
                    vi_new[i,:,k] = vi_y_new # Store values in new variable
        else:
            vi_new = var

        for i in range(nx):
            for j in range(ny_work):
                for k in range(nz):
                    gradphi_work = np.array((gradphi[i,j,k,0] , gradphi[i,j,k,1] , gradphi[i,j,k,2]))
                    vel_drift = np.cross(y_unit_3d[j],gradphi_work)/b_norm
                    vel_parallel = np.dot(y_unit_3d[j] , vi_new[i,j,k] ) #INSERT VPAR IN
                    velocity = vel_drift + vel_parallel
                    vectorx[i,j,k] = velocity[0]
                    vectory[i,j,k] = velocity[1]
                    vectorz[i,j,k] = velocity[2]
                    
        
        #Convert coordinate section        
        vrbl = visual.vtk_var(var_new,nx,ny_work,nz) # vtk variable
        vtk_path = visual.write_vtk_vector2(name,pts,vrbl,vectorx,vectory,vectorz,q)
        print "At t = %d, %d Steps remaining" % (q , ((t- q)/skip)) # Progress indicator
        q += skip
    
    #Write the Max and min values to file
#    mm_array = np.array((max,min))
#    np.savetxt('max_min_' + name + '.txt',mm_array)
    return
                

#==============================================================================
# Vector writing functions 
#==============================================================================

def cylinder(name, time, pi_fr = (2./3.), step = 0.5 , path = None, skip = 1):
    """
    Cylinder function writes the vector for a variable that is parallel along a field line (i.e Vpar)
    The variable is imported, interpolated, a mesh is created for the Cylindrical geometry
    unit vectors along the field direction are created using the y_unit_vector function then using imported variable
    the vectors are created and plotted
    Output: time stamped vtk files that contain the scalar variable and a vector variable of the imported 
    
    Inputs
        name: variable name to import (str)
        time: end time slice to convert to (int), if -1 entire dataset is imported
        step: The gap between the y slices for the y interpolation (float < 1) if 0 then raw data displayed.
        path: File path to data, if blank then the current working directory is used to look for data
        skip: The gap between time slices to convert.
        
    Outputs
        Timestamped vtk files of the specified variable, scalar and vector forms.
        A text file with the maximum and minimum files.   
    """    
    #Get the working dir
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    # Get the grid file name
    grid_file = visual.get("BOUT.inp","grid")

    # Find dimensions of data
    max_t, var_0, nx, ny, nz = visual.dim_all(name)
    ny_work = ny
    # Import coordinates from gridfile
    r = visual.nc_var(grid_file,'Rxy') # toroidal coordinates
    z = visual.nc_var(grid_file,'Zxy') # toroidal coordinates
        
    #Interpolate setup
    #Get number of new points
    if step != 0:
        ny_work = len(np.arange(0,(ny-1),step))
        r2 = visual.intrp_grd(ny,r,ny_work,step)
        z2 = visual.intrp_grd(ny,z,ny_work,step)
        r = r2
        z = z2
    #Define empty array with increased number of y values
    var_new = np.empty((nx,ny_work,nz),dtype=float)
    pts = visual.cylinder(r, z, nx, ny_work, nz, pi_fr) # vtk grid points
    
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
    max = np.zeros(t , dtype = float)
    min = np.zeros(t , dtype = float)
    
    #Find the y_unit vectors
    unitx, unity, unitz = y_unit_vector(pts, nx , ny_work , nz)
    #Create the empty vector variables
    vectorx = np.zeros((nx,ny_work,nz),dtype = float)
    vectory = np.zeros((nx,ny_work,nz),dtype = float)
    vectorz = np.zeros((nx,ny_work,nz),dtype = float)
    
    var_all = visual.collect(name)    
    
    q = 0
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    while q <= t-1:
            var = var_all[q] # collect variable
    
            #Find the min and max values
            max[q] = np.amax(var)
            min[q] = np.amin(var)
            
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
                
            for i in range(nx):
                for j in range(ny_work):
                    for k in range(nz):
                        vectorx[i,j,k] = unitx[i,j,k] * var_new[i,j,k]
                        vectory[i,j,k] = unity[i,j,k] * var_new[i,j,k]
                        vectorz[i,j,k] = unitz[i,j,k] * var_new[i,j,k]
            #Convert coordinate section        
            vrbl = visual.vtk_var(var_new,nx,ny_work,nz) # vtk variable
            vtk_path = visual.write_vtk_vector2(name,pts,vrbl, vectorx, vectory, vectorz ,q) # write vtk file
            print "At t = %d, %d Steps remaining" % (q,((t- q)/skip)) # Progress indicator
            q+= skip
            
    #Write the Max and min values to file
    mm_array = np.array((np.amax(max) , np.amin(min)))
    np.savetxt('max_min_' + name + '.txt',mm_array)
    return





def torus(name, time, step = 0.5, path = None, skip = 1, R = None, r = None , dr = None , Bt = None, q = None , isttok = False, default = False):
    """
    Torus function writes the vector for a variable that is parallel along a field line (i.e Vpar)
    The variable is imported, interpolated, a mesh is created for the Torous geometry
    unit vectors along the field direction are created using the torus_y_unit function then using imported variable
    the vectors are created and plotted
    Output: time stamped vtk files that contain the scalar variable and a vector variable of the imported 
    
    Inputs
        name: variable name to import (str)
        time: end time slice to convert to (int), if -1 entire dataset is imported
        step: The gap between the y slices for the y interpolation (float < 1) if 0 then raw data displayed.
        path: File path to data, if blank then the current working directory is used to look for data
        skip: The gap between time slices to convert.
        
        Torus settings, if left blank defaults are used all are floats.
        R: Major radius
        r: Minor radius
        dr: Radial width of domain
        Bt: Toroidal magnetic field
        q: Safety factor
        
        isttok: Use the ISTTOK specifications, (Boolean) R, r , dr , Bt , q = 0.46, 0.085, 0.02, 0.5, 5
        default: USe the Default specifications, (Boolean) R, r , dr , Bt , q = 2.0,  0.2, 0.05, 1.0, 5.0
        
    Outputs
        Timestamped vtk files of the specified variable, scalar and vector forms.
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
        
    if default == True: # Use default values?
        R,r,dr,Bt,q = 2.0 , 0.2 , 0.05, 1.0, 5.0

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
    
    #Initial max and min values
    max = np.zeros(t , dtype = float)
    min = np.zeros(t , dtype = float)
    
    #Find the y_unit vectors
    y_unit = torus_y_unit(ny_work)
    
    var_all = visual.collect(name)
    
    q = 0 
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    while q <= t-1:
        var = var_all[q] # collect variable
#        var = np.zeros((nx,ny_work,nz))
#        for i in range(nx):
#            for j in range(ny):
#                for k in range(nz):
#                    var[i,j,k] = 1
        
        #Find the min and max values
        max[q] = np.amax(var)
        min[q] = np.amin(var)
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
                
        # Assign values
        for i in range(nx):
            for j in range(ny_work):
                for k in range(nz):
                    vectorx[i,j,k] = var_new[i,j,k] * y_unit[j][0]
                    vectory[i,j,k] = var_new[i,j,k] * y_unit[j][1]
                    vectorz[i,j,k] = var_new[i,j,k] * 0 # zDirection is orthogonal to field therefore no component in z
        
        #Convert coordinate section        
        vrbl = visual.vtk_var(var_new,nx,ny_work,nz) # vtk variable
        vtk_path = visual.write_vtk_vector2(name,pts,vrbl,vectorx,vectory,vectorz,q)
        print "At t = %d, %d Steps remaining" % (q , ((t- q)/skip)) # Progress indicator
        q += skip
    
    #Write the Max and min values to file
    mm_array = np.array((np.amax(max) , np.abs(np.amin(min))))
    np.savetxt('max_min_' + name + '.txt',mm_array)
    return

# the elm function imports

def elm(name , time ,zShf_int_p = 0.25, path = None, skip = 1):
    """
       
    elm function: 
    Imports a variable that is parallel along a field line (for example; jpar) creates a mesh for
    the elm geoemtry, calculates the unit vectors along the field lines from the mesh (as this accounts for zShift).
    The function will output the vector and scalar form of this variable.
    
    
    
    Inputs
        name: variable name imported from BOUT.dmp files
        time: The end time slice that will be converted, if -1 entire dataset is imported
        zShift Interpolation Percentage: This is used to set the tolerance values for the linear interpolation.
        path: File path to data, if blank then the current working directory is used to look for data
    
    Output
        Timestamped vtk files of the specified variable, with scalar and vector variables.
        A text file with the maximum and minimum files.
                
    """
    zShf_int_p = visual.zShf_p_check(zShf_int_p)

    #Get the working dir
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    grid_file = visual.get("BOUT.inp","grid") # Import grid file name 

    # Get the dimensions and initial value of the variable data.    
    max_t, var_0, nx, ny, nz = visual.dim_all(name)
  
    #ELM
    #Import r,z and zshift from grid file
    r = visual.nc_var(grid_file,'Rxy') # ELM Coordinates
    z = visual.nc_var(grid_file,'Zxy') # ELM Coordinates
    zshift = visual.nc_var(grid_file,'zShift') #ELM Coordinates
    
    max_z,min_z = visual.z_shift_mm(nx,ny,zshift) #Find the max and minimum zshift values
    z_tol = min_z + ((max_z - min_z)* zShf_int_p) # Set tolerance value
    
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
    
    # Create max and min arrays
    max = np.zeros(t, dtype = float)
    min = np.zeros(t, dtype = float)

    unitx, unity, unitz = elm_field_vector(pts2,nx,ny2,nz) # Create unit vector map 
    
    vectorx = np.zeros((nx,ny2,nz),dtype = float) # Create the empty vector variables
    vectory = np.zeros((nx,ny2,nz),dtype = float)
    vectorz = np.zeros((nx,ny2,nz),dtype = float)
     #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    q = 0 

    var_all = visual.collect(name)    
    
    while q <= t-1:
        var = var_all[q] # collect variable

        # Find min and max values
        max[q] = np.amax(var)
        min[q] = np.amin(var)
    
        var2,ny3 = visual.zshift_interp3d(nx,ny,nz,zshift,z_tol,var) # Interpolate Variable
        
        for i in range(nx): 
            for j in range(ny2):
                for k in range(nz):
                    vectorx[i,j,k] = unitx[i,j,k] * var2[i,j,k]
                    vectory[i,j,k] = unity[i,j,k] * var2[i,j,k]
                    vectorz[i,j,k] = unitz[i,j,k] * var2[i,j,k]
        
        
        #Convert coordinate section
        vrbl = visual.vtk_var(var2,nx,ny2,nz) #vtk variable
        vtk_path = visual.write_vtk_vector2(name,pts2,vrbl,vectorx,vectory,vectorz,q)
        print "At t = %d, %d Steps remaining" % (q,((t- q)/skip))# Progress indicator
        q += skip
    print ""
    
    #Write the Max and min values to file
#    mm_array = np.array((np.median(max),np.median(min))
    mm_array = np.array((np.amax(max) , np.abs(np.amin(min))))
    np.savetxt('max_min_' + name + '.txt',mm_array) 
    return


