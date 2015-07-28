#==============================================================================
# eigen.py 
# Contains the function to plot the eigenvals on a graph then on the user clicking
# a value the eigenvalues associated with that point are written to a file.
#==============================================================================


# import relevant libraries
from boutdata import collect
import matplotlib.pyplot as plt
from numpy import argmin, amax, amin, arange
import sys
import numpy as np

#Import visual.py and visit
import settings as set
sys.path.insert(0,set.visit_dir)
import visit

import visual

# Define a fuction that will plot the eigenvalues
def plot_eigenvals(eigs, nx, ny, nz ,zshift, z_tol, ny2, name, pts ,data=None):
  fig = plt.figure() # Create a figure

  if data is None:
    # If no data supplied, only plot eigenvalues
    ax = fig.add_subplot(111)
  else:
    # If data, plot two figures
    ax = fig.add_subplot(211)
    
    # Check that the data has the right size
    if len(data.shape) != 4:
      raise ValueError("Expecting data to be 2D")
    if data.shape[0] != len(eigs):
      raise ValueError("First dimension of data must match length of eigs")


  eigs_r = eigs[:-1:2] # Collect the real eigenvalues
  eigs_i = eigs[1::2]	# Collect the imaginary eignevalues
  
  range_r = amax(eigs_r) - amin(eigs_r)
  range_i = amax(eigs_i) - amin(eigs_i)
  
  ax.plot(eigs_r, eigs_i, 'x') #Plot the real and imaginary eigenvalues
  ax.set_xlabel("Real component")
  ax.set_ylabel("Imaginary component")
  
  overplot, = ax.plot([], [], 'ok') #Create an empty plot
  
  if data is not None: #If there is some data plot it
    # Add a data plot
    ax2 = fig.add_subplot(212)
    vector_r, = ax2.plot([],[], '-k', label="Real")
    vector_i, = ax2.plot([],[], '-r', label="Imag")
    ax2.set_xlabel("X")
    ax2.set_ylabel("Eigenvector")
    ax2.legend()
	
  #Unpack variables
#  nx,ny,nz,zshift,z_tol,ny2,name,pts = data_variables
  def onclick(event):
    # Check if user clicked inside the plot
    if event.xdata is None:
      return

    # Find closest data point, but stretch axes so 
    # real and imaginary components are weighted equally
    dist = ((eigs_r - event.xdata)/range_r)**2 + ((eigs_i - event.ydata)/range_i)**2
    
    ind = argmin(dist)
    
    # Update the highlight plot
    overplot.set_data([eigs_r[ind]], [eigs_i[ind]])
    
    print("Eigenvalue number: %d (%e,%e)" % (ind, eigs_r[ind], eigs_i[ind]))
    
    if data is not None:

      #Draw the vtk file
      #Set the variables up
      var_r = data[2*ind,:] 
      var_i = data[2*ind+1,:]
      nx,ny,nz = visual.dimd(name)
     #Interpolate the variables
      var_r2,ny3 = visual.zshift_interp3d(nx,ny,nz,zshift,z_tol,var_r)
      var_i2,ny3 = visual.zshift_interp3d(nx,ny,nz,zshift,z_tol,var_i)
      #create the vtk variables
      vrbl_r = visual.vtk_var(var_r2,nx,ny2,nz)
      vrbl_i = visual.vtk_var(var_i2,nx,ny2,nz)
      #Write the vtk file
      vtk_path = visual.write_vtk_2(name,pts,vrbl_r,vrbl_i,ind)
      print "Eigenvalue %d written to vtk" % ind

    fig.canvas.draw()
  
  cid = fig.canvas.mpl_connect('button_press_event', onclick)
  plt.show()
      
# Draw function calls the previous functions and collects the eigenvalue data,
# writes x,y graph of the eigenvectors and then on clicking an eigen value the
# eigenvalues are written to the vtk format.
def draw(name,zshift_tol = 0.25):    
    grid_file = str(raw_input('Enter Grid File : '))
    eigs = collect("t_array")
    data = collect(name)
    #import the variable, set the nx,ny,nz values
    #write the pts variable
    nx,ny,nz = visual.dimd(name)
    #Import the R,Z values from the grid file
    r = visual.nc_var(grid_file,'Rxy') # R
    z = visual.nc_var(grid_file,'Zxy') # Z
    #  print r.shape,z.shape
    zshift = visual.nc_var(grid_file, 'zShift')	 # zShift
    z_tol = visual.z_shift_tol(nx,ny,zshift,zshift_tol)
    #Interpolate r,z,zshift values
    r2, ny2 = visual.zshift_interp2d(nx,ny,zshift,zshift_tol,r)
    z2, ny2 = visual.zshift_interp2d(nx,ny,zshift,zshift_tol,z)
    zshift2, ny2 = visual.zshift_interp2d(nx,ny,zshift,zshift_tol,zshift)

    pts2 = visual.elm(r2,z2,zshift2,nx,ny2,nz)
#    Plot eigen values and wri
    plot_eigenvals(eigs, nx, ny, nz ,zshift, z_tol, ny2, name, pts2, data=data)
      
      
#        name = str(raw_input('Enter variable name : '))
#        grid_file = str(raw_input('Enter Grid File : '))
#        eigs = collect("t_array")
#        data = collect(name)
#        #import the variable, set the nx,ny,nz values
#        #write the pts variable
#        nx,ny,nz = visual.dimd(name)
#        #Import the R,Z values from the grid file
#        r = visual.nc_var(grid_file,'Rxy') # R
#        z = visual.nc_var(grid_file,'Zxy') # Z
#        #  print r.shape,z.shape
#        zshift = visual.nc_var(grid_file, 'zShift')	 # zShift
#        tol = 0.25 # Set the tolerance value (change this to be imported from settings?
#        z_tol = visual.z_shift_tol(nx,ny,zshift,tol)
#        #Interpolate r,z,zshift values
#        r2, ny2 = visual.zshift_interp2d(nx,ny,zshift,z_tol,r)
#        z2, ny2 = visual.zshift_interp2d(nx,ny,zshift,z_tol,z)
#        zshift2, ny2 = visual.zshift_interp2d(nx,ny,zshift,z_tol,zshift)
#    
#        pts2 = visual.elm(r2,z2,zshift2,nx,ny2,nz)
#        #Package variables
#        #  data_variables = nx,ny,nz,zshift,z_tol,ny2,name,pts
#        plot_eigenvals(eigs, nx, ny, nz ,zshift, z_tol, ny2, name, pts2, data=data)