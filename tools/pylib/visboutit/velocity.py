# -*- coding: utf-8 -*-
"""
velocity.py
This library contains functions for finding; unit vectors, grad of a field, cross product, dot product, vector sum and write velocity vectors in the cylinder, torus and elm coorinate systems.
"""
import numpy as np
import visual
from scipy import interpolate
import os
import circle

# This Function calulcates y unit vectors from a mesh
def y_vector_unit(pts,nx,ny,nz):
    xcoord, ycoord, zcoord = pts
    # Create the unit vector variable
    unit_vector_x = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_y = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_z = np.zeros((nx,ny,nz),dtype = float)
    
    # For every x,y,z point calculate a unit vector along the field line
    # by calculating the vector between two points from mesh (this already accounts for zShift)
    print 'Creating y vector map:'
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

# This Function calulcates x unit vectors from a mesh
def x_vector_unit(pts,nx,ny,nz):
    xcoord, ycoord, zcoord = pts
    # Create the unit vector variable
    unit_vector_x = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_y = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_z = np.zeros((nx,ny,nz),dtype = float)
    
    # For every x,y,z point calculate a unit vector along the field line
    # by calculating the vector between two points from mesh (this already accounts for zShift)
    print 'Creating x vector map:'
    for i in range(nx-1):
        for j in range(ny-1):
            percent = (float(i)/(float(nx - 2)))*100
            print '%d %% Complete' % percent, "   \r", # Progress percent
            for k in range(nz-1):
                r1 = np.array((xcoord[ i , j , k ] , ycoord[ i , j , k ] , zcoord[ i , j , k ])) # Vector of initial point
                r2 = np.array((xcoord[ i+1 , j , k ] , ycoord[ i+1 , j , k ]  , zcoord[ i+1 , j , k ])) # Vector of secondary point
                r = r2-r1
                unit = r/ np.linalg.norm(r) # Normalise vector
                # Assign the vectors
                unit_vector_x[i,j,k] = unit[0]
                unit_vector_y[i,j,k] = unit[1]
                unit_vector_z[i,j,k] = unit[2]
    print "" # For closing updated print statement
    return unit_vector_x, unit_vector_y, unit_vector_z

# This Function calulcates z unit vectors from a mesh
def z_vector_unit(pts,nx,ny,nz):
    xcoord, ycoord, zcoord = pts
    # Create the unit vector variable
    unit_vector_x = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_y = np.zeros((nx,ny,nz),dtype = float)
    unit_vector_z = np.zeros((nx,ny,nz),dtype = float)
    
    # For every x,y,z point calculate a unit vector along the field line
    # by calculating the vector between two points from mesh (this already accounts for zShift)
    print 'Creating z vector map:'
    for i in range(nx-1):
        for j in range(ny-1):
            percent = (float(i)/(float(nx - 2)))*100
            print '%d %% Complete' % percent, "   \r", # Progress percent
            for k in range(nz-1):
                r1 = np.array((xcoord[ i , j , k ] , ycoord[ i , j , k ] , zcoord[ i , j , k ])) # Vector of initial point
                r2 = np.array((xcoord[ i , j , k+1 ] , ycoord[ i , j , k+1 ]  , zcoord[ i , j , k+1 ])) # Vector of secondary point
                r = r2-r1
                unit = r/ np.linalg.norm(r) # Normalise vector
                # Assign the vectors
                unit_vector_x[i,j,k] = unit[0]
                unit_vector_y[i,j,k] = unit[1]
                unit_vector_z[i,j,k] = unit[2]
    print "" # For closing updated print statement
    return unit_vector_x, unit_vector_y, unit_vector_z

# This function calculates grad of a variable in a particular mesh
def grad_func(pts,var,nx,ny,nz):
    xcoord, ycoord, zcoord = pts
    
    grad_i = np.zeros((nx, ny, nz), dtype=float)
    grad_j = np.zeros((nx, ny, nz), dtype=float)
    grad_k = np.zeros((nx, ny, nz), dtype=float)

    print 'Computing Grad:'
    for i in range(nx-2):
        for j in range(ny-2):
            percent = (float(i)/(float(nx - 3)))*100
            print '%d %% Complete' % percent, "   \r", # Progress percent
            for k in range(nz-2):
                x,y,z = i+1 , j+1 , k+1
                #Import values
                r0 = np.array((xcoord[x, y, z], ycoord[x, y, z], zcoord[x, y, z]), dtype=float)
                r1 = np.array((xcoord[x+1, y, z], ycoord[x+1, y, z], zcoord[x+1, y, z]), dtype=float)
                r2 = np.array((xcoord[x, y+1, z], ycoord[x, y+1, z], zcoord[x, y+1, z]), dtype=float)
                r3 = np.array((xcoord[x, y, z+1], ycoord[x, y, z+1], zcoord[x, y, z+1]), dtype=float)
                r4 = np.array((xcoord[x-1, y, z], ycoord[x-1, y, z], zcoord[x-1, y, z]), dtype=float)
                r5 = np.array((xcoord[x, y-1, z], ycoord[x, y-1, z], zcoord[x, y-1, z]), dtype=float)
                r6 = np.array((xcoord[x, y, z-1], ycoord[x, y, z-1], zcoord[x, y, z-1]), dtype=float)
                phi0 = var[x, y, z]
                phi1 = var[x+1, y, z]
                phi2 = var[x, y+1, z]
                phi3 = var[x, y, z+1]
                phi4 = var[x-1, y, z]
                phi5 = var[x, y-1, z]
                phi6 = var[x, y, z-1]
                
                # Diff in points
                x_displ_pve = r1 - r0
                y_displ_pve = r2 - r0
                z_displ_pve = r3 - r0 
                x_displ_nve = r4 - r0
                y_displ_nve = r5 - r0
                z_displ_nve = r6 - r0 
                
                # solve for dphi/dx, dphi/dy, dphi/dz in +ve direction
                mtrx_pve_cnsts = np.array([[1 , 0 , 0 , 0 ] ,
                                           [1 , (x_displ_pve[0]) , (x_displ_pve[1]) , (x_displ_pve[2]) ],
                                            [1, (y_displ_pve[0]) , (y_displ_pve[1]) , (y_displ_pve[2]) ],
                                            [1, (z_displ_pve[0]) , (z_displ_pve[1]) , (z_displ_pve[2]) ]])                
                mtrx_nve_cnsts = np.array([[1 , 0 , 0 , 0 ] ,
                                           [1 , (x_displ_nve[0] ) , (x_displ_nve[1]) , (x_displ_nve[2]) ],
                                            [1, (y_displ_nve[0]) , (y_displ_nve[1]) , (y_displ_nve[2]) ],
                                            [1, (z_displ_nve[0]) , (z_displ_nve[1]) , (z_displ_nve[2]) ]])  
                mtrx_phi_pve = np.array([[phi0] , [phi1] , [phi2] , [phi3]])
                mtrx_phi_nve = np.array([[phi0] , [phi4] , [phi5] , [phi6]])
                mtrx_soln_pve = np.linalg.solve(mtrx_pve_cnsts , mtrx_phi_pve)
                mtrx_soln_nve = np.linalg.solve(mtrx_nve_cnsts , mtrx_phi_nve)
                
                dphi_dx = np.average((mtrx_soln_pve[1],mtrx_soln_nve[1]))
                dphi_dy = np.average((mtrx_soln_pve[2],mtrx_soln_nve[2]))
                dphi_dz = np.average((mtrx_soln_pve[3],mtrx_soln_nve[3]))
                

                grad_i[x,y,z] = dphi_dx
                grad_j[x,y,z] = dphi_dy
                grad_k[x,y,z] = dphi_dz
    print ""
    return grad_i, grad_j, grad_k

# This function calculates the vector cross product between two vectors
def cross_product(avector, bvector, nx, ny, nz):
    print "Calculating cross product"
    avector_i, avector_j, avector_k = avector
    bvector_i, bvector_j, bvector_k = bvector
    cp_i = np.zeros((nx, ny, nz), dtype=float)
    cp_j = np.zeros((nx, ny, nz), dtype=float)
    cp_k = np.zeros((nx, ny, nz), dtype=float)
    for i in range (nx):
        for j in range(ny):
            percent = (float(i)/float(nx-1)) * 100
            print '%d %% Complete' % percent, "   \r", # Progress percent
            for k in range(nz-1):
                avector_work = np.array((avector_i[i,j,k], avector_j[i,j,k], avector_k[i,j,k]))
                bvector_work = np.array((bvector_i[i,j,k], bvector_j[i,j,k], bvector_k[i,j,k]))
                cp_work = np.cross(avector_work, bvector_work)
                cp_i[i,j,k] = cp_work[0]
                cp_j[i,j,k] = cp_work[1]
                cp_k[i,j,k] = cp_work[2]
    print ""
    return cp_i, cp_j, cp_k

# This function calculates the dot product of a scalar and a vector
def dot_product(ascalar, bvector, nx, ny, nz):
    print "Calculating dot product"
    bvector_i, bvector_j, bvector_k = bvector
    dp_i = np.zeros((nx, ny, nz), dtype=float)
    dp_j = np.zeros((nx, ny, nz), dtype=float)
    dp_k = np.zeros((nx, ny, nz), dtype=float)
    for i in range (nx):
        for j in range(ny):
            percent = (float(i)/float(nx-1)) * 100
            print '%d %% Complete' % percent, "   \r", # Progress percent
            for k in range(nz-1):
                ascalar_work = np.array(ascalar[i,j,k])
                bvector_work = np.array((bvector_i[i,j,k], bvector_j[i,j,k], bvector_k[i,j,k]))
                dp_work = np.dot(ascalar_work, bvector_work)
                dp_i[i,j,k] = dp_work[0]
                dp_j[i,j,k] = dp_work[1]
                dp_k[i,j,k] = dp_work[2]
    print ""
    return dp_i, dp_j, dp_k

# This function calculates the sum of two vectors
def vec_add(avector, bvector, nx, ny, nz):
    print "Calculating vec sum"
    avector_i, avector_j, avector_k = avector
    bvector_i, bvector_j, bvector_k = bvector
    add_i = np.zeros((nx, ny, nz), dtype=float)
    add_j = np.zeros((nx, ny, nz), dtype=float)
    add_k = np.zeros((nx, ny, nz), dtype=float)
    for i in range (nx):
        for j in range(ny):
            percent = (float(i)/float(nx-1)) * 100
            print '%d %% Complete' % percent, "   \r", # Progress percent
            for k in range(nz-1):
                avector_work = np.array((avector_i[i,j,k], avector_j[i,j,k], avector_k[i,j,k]))
                bvector_work = np.array((bvector_i[i,j,k], bvector_j[i,j,k], bvector_k[i,j,k]))
                add_work = avector_work + bvector_work
                add_i[i,j,k] = add_work[0]
                add_j[i,j,k] = add_work[1]
                add_k[i,j,k] = add_work[2]
    print ""
    return add_i, add_j, add_k

def velvar_check(var_list):
    if 'phi' in var_list:
        phi_name = 'phi'
    else:
        phi_name = str(raw_input('Name of potential variable (i.e. phi): '))
        
    if 'Tnorm' in var_list:
        phi_norm_name = 'Tnorm'
    else:
        phi_norm_name = str(raw_input('Name of potential normalisation variable (i.e.Tnorm): ' ))
 
    if 'Bnorm' in var_list:
        B_norm_name = 'Bnorm'
    else:
        B_norm_name = str(raw_input('Name of Magnetic normalisation variable (i.e. Bnorm): '))
         
    if 'Vi' in var_list:
        vi_name = 'Vi'
        NVi_name = None
        Ne_name = None
    else:
        vi_name = None
        if 'NVi' in var_list:
            NVi_name = 'NVi'
        else:
            NVi_name = str(raw_input('Name of NVi variable: '))
        if 'Ne' in var_list:
            Ne_name = 'Ne'
        else:
            Ne_name = str(raw_input('Name of Ne variable: '))
             
    if 'Cs0' in var_list:
        Cs0_name = 'Cs0'
    else:
        Cs0_name = str(raw_input('Name of normalised sound speed (i.e. Cs0): '))

    return phi_name, phi_norm_name, B_norm_name, vi_name, NVi_name, Ne_name, Cs0_name
    

# This function calculates and writes velocity vectors to the VTK format for the cylinder coordinates
def cylinder(time = -1, pi_fr = (2./3.), step = 0.5 , path = None, skip = 1):
    
    # Variable setup    
    var_list = visual.var_list()
    
    phi_name, phi_norm_name, B_norm_name, vi_name, NVi_name, Ne_name, Cs0_name = velvar_check(var_list)
    
    #Get the working dir
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    # Get the grid file name
    grid_file = visual.get("BOUT.inp","grid")

    # Find dimensions of data
    max_t, var_0, nx, ny, nz = visual.dim_all(phi_name)
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
    phi_new = np.empty((nx,ny_work,nz),dtype=float)
    vi_new = np.empty((nx,ny_work,nz),dtype=float)
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
    
#    #Initial max and min values
#    max = np.zeros(t , dtype = float)
#    min = np.zeros(t , dtype = float)
    
    #Find the unit vectors
#    x_vector_i, x_vector_j, x_vector_k = x_vector_unit(pts, nx, ny_work, nz)
#    x_vector = x_vector_i, x_vector_j, x_vector_k
    
    y_vector_i, y_vector_j, y_vector_k = y_vector_unit(pts, nx, ny_work, nz)
    y_vector = y_vector_i, y_vector_j, y_vector_k
    
#    z_vector_i, z_vector_j, z_vector_k = z_vector_unit(pts, nx, ny_work, nz)
#    z_vector = z_vector_i, z_vector_j, z_vector_k
    
    # Import normalisation factors
    B_norm = visual.collect(B_norm_name)
    phi_norm = visual.collect(phi_norm_name)
    cs0 = visual.collect(Cs0_name)

    
    # Import data
    phi_all = visual.collect(phi_name) * phi_norm #phi data
    if vi_name == None:
        ne = visual.collect(Ne_name)
        nvi = visual.collect(NVi_name)
        vi_all = (nvi/ne) *cs0 #vi data        
    else:
        vi_all = visual.collect(vi_name) * cs0
    
    #Create the empty vector variables
    grad_phix = np.zeros((nx,ny_work,nz),dtype=float)
    grad_phiy, grad_phiz = grad_phix, grad_phix    
    
    q = 0
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    
    while q <= t-1:
        phi = phi_all[q]
        vi = vi_all[q]
        if ny != ny_work: #Interpolate phi for all y values
            for i in range(nx):
                for k in range(nz):
                    phi_y = phi[i,:,k]
                    y = np.arange(0,ny) # Arrange y values in array
                    f = interpolate.interp1d(y,phi_y) # Interpolate the data
                    y_new = np.arange(0,(ny-1),step) # Arrange new y values 
                    var_y_new = f(y_new) # interpolate y values
                    phi_new[i,:,k] = var_y_new # Store values in new variable
        else:
            phi_new = phi
        
        if ny != ny_work: #Interpolate vi for all y values
            for i in range(nx):
                for k in range(nz):
                    vi_y = vi[i,:,k]
                    y = np.arange(0,ny) # Arrange y values in array
                    f = interpolate.interp1d(y,vi_y) # Interpolate the data
                    y_new = np.arange(0,(ny-1),step) # Arrange new y values 
                    var_y_new = f(y_new) # interpolate y values
                    vi_new[i,:,k] = var_y_new # Store values in new variable
        else:
            vi_new = vi

        # Calculate grad_phi
        grad_phix, grad_phiy, grad_phiz = grad_func(pts, phi_new, nx, ny_work, nz)
        grad_phi = grad_phix, grad_phiy, grad_phiz
        
        #Setup arrays for magnetic contribution to velocity
        mag_vx = np.zeros((nx, ny_work, nz), dtype = float)        
        mag_vy, mag_vz = mag_vx, mag_vx
        
        #Calculate magnetic contribution to velocity.
        mag_vx, mag_vy, mag_vz = cross_product(y_vector, grad_phi, nx, ny_work, nz)
        
        # Normalise magnetic vel
        mag_vx, mag_vy, mag_vz = (mag_vx/B_norm), (mag_vy/B_norm), (mag_vz/B_norm)
        mag_v = mag_vx, mag_vy, mag_vz
       
       # Setup arrays for velocity
        velx = np.zeros((nx, ny_work, nz), dtype = float)
        vely, velz = velx, velx
        
        # Calcuate Vi parallel vector
        elec_vel = dot_product(vi_new, y_vector, nx, ny_work, nz)
        
        # Calculate velocity
        velx, vely, velz = vec_add(elec_vel, mag_v, nx, ny_work, nz)

        vrbl = visual.vtk_var(phi_new,nx,ny_work,nz) # vtk variable

        vtk_path = visual.write_vtk_vector2('velocity', pts, vrbl, velx, vely, velz, q)
        print "At t = %d, %d Steps remaining" % (q,((t- q)/skip)) # Progress indicator
        q += skip
    return
    
# This function calculates and writes velocity vectors to the VTK format for the torus coordinates
def torus(time = -1, step = 0.5, skip = 1 , path = None, R = None, r = None , dr = None , Bt = None, q = None , isttok = True , default = False):
  
    # Variable setup    
    var_list = visual.var_list()
    
    phi_name, phi_norm_name, B_norm_name, vi_name, NVi_name, Ne_name, Cs0_name = velvar_check(var_list)
    

    #Get the working dir
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    # Get the grid file name
    grid_file = visual.get("BOUT.inp","grid")

    # Find the max_t, initial value and shape of variable
    max_t, var_0, nx, ny, nz = visual.dim_all(phi_name)
    ny_work = ny
    #Interpolate setup
    #Get number of new points
    if step != 0:
        ny_work = len(np.arange(0,(ny-1),step))
    
    if isttok == True: # Use ISTTOK specifications?
        R,r,dr,Bt,q = 0.46, 0.085, 0.02, 0.5, 5
        
    if default == True: # Use default values?
        R,r,dr,Bt,q = 2.0 , 0.2 , 0.05, 1.0, 5.0
        
    # Import coordinates from function in circle library
    r,z = circle.coordinates(nx,ny_work,nz,R,r,dr,Bt,q) # toroidal coordinates
    pts = visual.torus(r, z, nx, ny_work, nz)  # vtk grid points
        
    #Interpolate setup
    #Get number of new points
    if step != 0:
        ny_work = len(np.arange(0,(ny-1),step))
   
    #Define empty array with increased number of y values
    phi_new = np.empty((nx,ny_work,nz),dtype=float)
    vi_new = np.empty((nx,ny_work,nz),dtype=float)
    
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

    #Find the unit vectors
#    x_vector_i, x_vector_j, x_vector_k = x_vector_unit(pts, nx, ny_work, nz)
#    x_vector = x_vector_i, x_vector_j, x_vector_k
    
    y_vector_i, y_vector_j, y_vector_k = y_vector_unit(pts, nx, ny_work, nz)
    y_vector = y_vector_i, y_vector_j, y_vector_k
    
#    z_vector_i, z_vector_j, z_vector_k = z_vector_unit(pts, nx, ny_work, nz)
#    z_vector = z_vector_i, z_vector_j, z_vector_k
    
    # Import normalisation factors
    phi_norm = visual.collect(phi_norm_name) 
    B_norm = visual.collect(B_norm_name)
    cs0 = visual.collect(Cs0_name)

    
    # Import data
    phi_all = visual.collect(phi_name) * phi_norm #phi data
    if vi_name != None:
        vi_all = visual.collect(vi_name) * cs0
    else:
        ne = visual.collect(Ne_name)
        nvi = visual.collect(NVi_name)
        vi_all = (nvi/ne) *cs0 #vi data
    
    #Create the empty vector variables
    grad_phix = np.zeros((nx,ny_work,nz),dtype=float)
    grad_phiy, grad_phiz = grad_phix, grad_phix    
    
    q = 0
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    
    while q <= t-1:
        phi = phi_all[q]
        vi = vi_all[q]
        if ny != ny_work: #Interpolate phi for all y values
            for i in range(nx):
                for k in range(nz):
                    phi_y = phi[i,:,k]
                    y = np.arange(0,ny) # Arrange y values in array
                    f = interpolate.interp1d(y,phi_y) # Interpolate the data
                    y_new = np.arange(0,(ny-1),step) # Arrange new y values 
                    var_y_new = f(y_new) # interpolate y values
                    phi_new[i,:,k] = var_y_new # Store values in new variable
        else:
            phi_new = phi
        
        if ny != ny_work: #Interpolate vi for all y values
            for i in range(nx):
                for k in range(nz):
                    vi_y = vi[i,:,k]
                    y = np.arange(0,ny) # Arrange y values in array
                    f = interpolate.interp1d(y,vi_y) # Interpolate the data
                    y_new = np.arange(0,(ny-1),step) # Arrange new y values 
                    var_y_new = f(y_new) # interpolate y values
                    vi_new[i,:,k] = var_y_new # Store values in new variable
        else:
            vi_new = vi

        # Calculate grad_phi
        grad_phix, grad_phiy, grad_phiz = grad_func(pts, phi_new, nx, ny_work, nz)
        grad_phi = grad_phix, grad_phiy, grad_phiz
        
        #Setup arrays for magnetic contribution to velocity
        mag_vx = np.zeros((nx, ny_work, nz), dtype = float)        
        mag_vy, mag_vz = mag_vx, mag_vx
        
        #Calculate magnetic contribution to velocity.
        mag_vx, mag_vy, mag_vz = cross_product(y_vector, grad_phi, nx, ny_work, nz)
        
        # Normalise magnetic vel
        mag_vx, mag_vy, mag_vz = (mag_vx/B_norm), (mag_vy/B_norm), (mag_vz/B_norm)
        mag_v = mag_vx, mag_vy, mag_vz
        
        # Setup arrays for velocity
        velx = np.zeros((nx, ny_work, nz), dtype = float)
        vely, velz = velx, velx
        
        # Calculate Vi parallel vectors
        elec_vel = dot_product(vi_new, y_vector, nx, ny_work, nz)

        # Calculate velocity
        velx, vely, velz = vec_add(elec_vel, mag_v, nx, ny_work, nz)

        vrbl = visual.vtk_var(phi_new,nx,ny_work,nz) # vtk variable
        
        vtk_path = visual.write_vtk_vector2('velocity', pts, vrbl, velx, vely, velz  , q)
        print "At t = %d, %d Steps remaining" % (q,((t- q)/skip)) # Progress indicator
        q += skip
    return

# This function calculates and writes velocity vectors to the VTK format for the cylinder coordinates
def elm(time = -1, zShf_int_p = 0.25, path = None, skip = 1):
   
    # Variable setup    
    var_list = visual.var_list()
    
    phi_name, phi_norm_name, B_norm_name, vi_name, NVi_name, Ne_name, Cs0_name = velvar_check(var_list)
    
    zShf_int_p = visual.zShf_p_check(zShf_int_p)
    #Get the working dir
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    # Get the grid file name
    grid_file = visual.get("BOUT.inp","grid")

    # Find the max_t, initial value and shape of variable
    max_t, var_0, nx, ny, nz = visual.dim_all(phi_name)

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

    #Define empty array with increased number of y values
    phi_new = np.empty((nx,ny2,nz),dtype=float)
    vi_new = np.empty((nx,ny2,nz),dtype=float)
    
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
    
#    #Initial max and min values
#    max = np.zeros(t , dtype = float)
#    min = np.zeros(t , dtype = float)
    
    #Find the unit vectors
#    x_vector_i, x_vector_j, x_vector_k = x_vector_unit(pts2, nx, ny2, nz)
#    x_vector = x_vector_i, x_vector_j, x_vector_k
    
    y_vector_i, y_vector_j, y_vector_k = y_vector_unit(pts2, nx, ny2, nz)
    y_vector = y_vector_i, y_vector_j, y_vector_k
    
#    z_vector_i, z_vector_j, z_vector_k = z_vector_unit(pts2, nx, ny2, nz)
#    z_vector = z_vector_i, z_vector_j, z_vector_k
    
    
    # Import normalisation factors
    phi_norm = visual.collect(phi_norm_name) * 1.602177e-19
    B_norm = visual.collect(B_norm_name)
    cs0 = visual.collect(Cs0_name)

    # Import data
    phi_all = visual.collect(phi_name) * phi_norm #phi data
    if vi_name != None:
        vi_all = visual.collect(vi_name) * cs0
    else:
        ne = visual.collect(Ne_name)
        nvi = visual.collect(NVi_name)
        vi_all = (nvi/ne) *cs0 #vi data
    
    #Create the empty vector variables
    velx = np.zeros((nx,ny2,nz),dtype = float)
    vely = np.zeros((nx,ny2,nz),dtype = float)
    velz= np.zeros((nx,ny2,nz),dtype = float)
    
    grad_phix = np.zeros((nx,ny2,nz),dtype=float)
    grad_phiy, grad_phiz = grad_phix, grad_phix    
    
    q = 0
    #For the entire t range import the spacial values, find min max values, interpolate y, coordinate transform, write to vtk
    
    while q <= t-1:
        phi = phi_all[q]
        vi = vi_all[q]
        
        phi_new,ny3 = visual.zshift_interp3d(nx,ny,nz,zshift,z_tol,phi) # Interpolate Variable
        vi_new,ny3 = visual.zshift_interp3d(nx,ny,nz,zshift,z_tol,vi) # Interpolate Variable

        # Calculate grad_phi
        grad_phix, grad_phiy, grad_phiz = grad_func(pts2, phi_new, nx, ny2, nz)
        grad_phi = grad_phix, grad_phiy, grad_phiz
        
        #Setup arrays for magnetic contribution to velocity
        mag_vx = np.zeros((nx, ny2, nz), dtype = float)        
        mag_vy, mag_vz = mag_vx, mag_vx
        
        #Calculate magnetic contribution to velocity.
        mag_vx, mag_vy, mag_vz = cross_product(y_vector, grad_phi, nx, ny2, nz)
        
        # Normalise magnetic vel
        mag_vx, mag_vy, mag_vz = (mag_vx/B_norm), (mag_vy/B_norm), (mag_vz/B_norm)
        mag_v = mag_vx, mag_vy, mag_vz
        # Setup arrays for velocity
        velx = np.zeros((nx, ny2, nz), dtype = float)
        vely, velz = velx, velx
        
        elec_vel = dot_product(vi_new, y_vector, nx, ny2, nz)
        
        # Calculate velocity
        velx, vely, velz = vec_add(elec_vel, mag_v, nx, ny2, nz)

        vrbl = visual.vtk_var(phi_new,nx,ny2,nz) # vtk variable

        vtk_path = visual.write_vtk_vector2('velocity', pts2, vrbl, velx, vely, velz, q)
        print "At t = %d, %d Steps remaining" % (q,((t- q)/skip)) # Progress indicator
        q += skip
    return

