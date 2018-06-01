"""
 eigen.py 
 Contains the function to plot the eigenvals on a graph then on the user clicking
 a value the eigenvalues associated with that point are written to a file.
"""

# import relevant libraries
from . import visual, vtk
from boutdata import collect
import matplotlib.pyplot as plt
from numpy import argmin, amax, amin
import os

# Define a fuction that will plot the eigenvalues
def plot_eigenvals(eigs, name, coord, zShf_int_p, t_step, R,r,dr,Bt,q, c_step, data=None):
    """
    This function is used to plot the eigenvectors and convert the specific sets of data. This function is called by the draw function and should not be called directly
    """
    fig = plt.figure() # Create a figure
    plt.title("A graph of imported BOUT++ eigenvectors")
    ax = fig.add_subplot(111) # Add a plot to the figure
    
    # Check that the data has the right size
    if len(data.shape) != 4:
        raise ValueError("Expecting data to be 2D")
    if data.shape[0] != len(eigs):
        raise ValueError("First dimension of data must match length of eigs")
        
    eigs_r = eigs[:-1:2] # Collect the real eigenvalues
    eigs_i = eigs[1::2]    # Collect the imaginary eignevalues
    
    range_r = amax(eigs_r) - amin(eigs_r)
    range_i = amax(eigs_i) - amin(eigs_i)
    
    ax.plot(eigs_r, eigs_i, 'x') #Plot the real and imaginary eigenvalues
    ax.set_xlabel("Real component")
    ax.set_ylabel("Imaginary component")    
    
    overplot, = ax.plot([], [], 'ok') #Create an empty plot

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
            #       Draw the vtk file
            #    Set the variables up
            var_r = data[2*ind,:]
            var_i = data[2*ind+1,:]

            if coord.lower() == 'elm': # Add in coordinate options
                vrbl_r , pts = vtk.elm_t(name , 2*ind , zShf_int_p)
                vrbl_i , pts = vtk.elm_t(name , 2*ind + 1, zShf_int_p)

            if coord.lower() == 'torus': # Toroidal coordinate system
                vrbl_r , pts = vtk.torus_t(name, 2*ind, t_step, R=R, r=r, dr=dr, Bt = Bt, q = q)
                vrbl_i , pts = vtk.torus_t(name, 2*ind+1, t_step, R=R, r=r, dr=dr, Bt = Bt, q = q)
                
            if coord.lower() == 'cylinder': # Cylindrical coordinate system
                vrbl_r , pts = vtk.cylinder_t(name, 2*ind, c_step)
                vrbl_i , pts = vtk.cylinder_t(name, 2*ind+1, c_step) 

            
            vtk_path = visual.write_vtk_2(name,pts,vrbl_r,vrbl_i,ind)
            print("Eigenvalue %d written to vtk" % ind)
        fig.canvas.draw()
        
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
      
# Draw function calls the previous functions and collects the eigenvalue data,
# writes x,y graph of the eigenvectors and then on clicking an eigen value the
# eigenvalues are written to the vtk format.
def draw(name, path = None):
    
    """
    The draw function handels the importing of the data, and setting up of the conversion before calling the plot_eigenvals to convert the data
    
    Inputs:
        name: variable name to be converted
        path: bath to raw BOUT++ files, if they are not within the current working directory
        
    Output:
        This function (which utilises the plot_eigenvals function) outputs converted eigenvalue data.
    
    
    """
    # Set the work_dir
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()

    # Create batch dir
    if not os.path.exists('batch'):
        os.makedirs("batch")
    
    # Create Image dir
    if not os.path.exists('images'):
        os.makedirs('images')

    eigs = collect("t_array") # Import Eigen values
    var = collect(name) # Import variable
    
    
    coord_check = False
    coord_list = ['cylinder','torus','elm']
    
    # Check if the coordinate system is one of the accepted systems.
    while coord_check == False:
        coord = str(raw_input('Enter coordinate system; cylinder, torus or elm: ' ))
        for i in coord_list:
            if coord.lower() == i.lower():
                coord_check = True
    
    zShf_int_p, t_step, R,r,dr,Bt,q, c_step = None, None, None, None, None, None, None, None

    # Prompt for the coordinate
    if coord.lower() == 'elm':
        zShf_int_p = float(raw_input('Enter zShift Interpolation Percent (Default 0.25) :'))
        
    if coord.lower() == 'torus':
        t_step = float(raw_input('Enter step value for Interpolation (Default 0.5): '))
        use_defaults = str(raw_input('Use default specifications (y/n)? '))
        if use_defaults.lower() == 'y':
            R,r,dr,Bt,q = 2.0, 0.2, 0.05, 1.0, 5.0
        else:
            use_isttok = str(raw_input('Use ISTTOK specifications R,r,dr,Bt,q = 0.46, 0.085, 0.02, 0.5, 5 (y/n):  '))
            if use_isttok.lower() == 'y':
                R,r,dr,Bt,q = 0.46, 0.085, 0.02, 0.5, 5
            else:
                print("Enter values for specifications")
                R = float(raw_input('R = '))
                r = float(raw_input('r = '))
                dr = float(raw_input('dr = '))
                Bt = float(raw_input('Bt = '))
                q = float(raw_input('q = '))
            
    if coord.lower() == 'cylinder':
        c_step = float(raw_input('Enter step value for Interpolation (Default 0.5): '))


#    # Plot eigen values and write the associated vtk file
    plot_eigenvals(eigs , name , coord, zShf_int_p, t_step, R,r,dr,Bt,q, c_step ,var)
