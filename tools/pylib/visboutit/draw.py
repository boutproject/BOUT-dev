"""
draw.py 
This library contains functions that display converted BOUT++ data in a VisIt window allowing the user to modify the plot before rendering out an image sequence of the plot.

There are three main functions that handle this; draw (scalar data), vector (field parallel vector data) and eigen (data from eigensolver). 


"""
from . import visual
import sys
import os
import numpy as np

# Import the VisIt library
sys.path.insert(0,visual.visit_dir)
import visit

#==============================================================================
# view_vector draws the vector variable withina VisIt window allowing the user 
# to setup the desired view settings.
#==============================================================================
def view_vector(work_dir,name,max,min):
    vtk_path = work_dir + "/batch/" + name + "_vector_*.vts database" # Set vtkfile path
    visit.OpenDatabase(vtk_path) # Open database
    visit.AddPlot("Vector",name + "_vector", 1, 1) #Draw a Pseudocolor Plot of the variable

#    VectorAtts = visit.VectorAttributes()
#    VectorAtts.nVectors = 4000 # Increase the number of vectors
#    visit.SetPlotOptions(VectorAtts)
    visit.DrawPlots() # Draw the Plots
    
    # Save the Visit Session
    session_name = raw_input('Enter a session file name:')
    session_path = work_dir + "/" + session_name + ".session"
    visit.SaveSession(session_path)
    # Close visit session
    visit.DeleteAllPlots()
    visit.CloseDatabase(vtk_path)
    return session_path,session_name

#==============================================================================
# Export an image sequence of the vector plot across the entire time range
#==============================================================================

def image_vector(session_path,img_dir,name,t,session_name,max_imp,min_imp,skip):
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
    sys.path.insert(0,visual.visit_dir)
    import visit
    # Load session and initialise at time 0
    visit.RestoreSession(session_path,0)
    visit.SetTimeSliderState(0)
    # Export an Image sequence of the variable for every time base
    for i in indicies:
        time = i * skip
        visit.SetTimeSliderState(i) # Change timer slider
        VectorAtts = visit.VectorAttributes()
#        VectorAtts.nVectors = 4000 # Increase the number of vectors
#        visit.SetPlotOptions(VectorAtts)
        # If user would like fixed max and mind then assign max and min values
        if max_imp != 0:
            VectorAtts.max = max_imp
            VectorAtts.maxFlag = 1
        if max_imp == 0:
            VectorAtts.maxFlag = 0
        if min_imp != 0:
            VectorAtts.min = min_imp
            VectorAtts.minFlag = 1
        if min_imp == 0:
            VectorAtts.minFlag = 0
        visit.SetPlotOptions(VectorAtts)
        visit.DrawPlots() # Draw plot
        # Save a png of the plot
        s = visit.SaveWindowAttributes()
        s.outputToCurrentDirectory = 0
        s.outputDirectory = outputdir
        s.family = 0
        s.fileName = '%s_%s_image_%04d' % (name,session_name,time)
        s.format = s.PNG
        s.width = visual.img_width
        s.height = visual.img_height
        visit.SetSaveWindowAttributes(s)
        visit.SaveWindow()
    # Close visit session
    visit.DeleteAllPlots()
    visit.Close()
    

#==============================================================================
# User Functions
#==============================================================================

# Draw the scalar variable and then save a session
def draw(name , skip = 1 , path = None):
    
    """
    This function imports a scalar variable from converted BOUT++ data and displays a Pseudocolor plot in a VisIt session
    
    Inputs:
        name: name of variable converted
        skip: skip value used in conversion
        path: path to data if not located within working directory.
        
    Ouput:
        Image sequence of the VisIt plot.
        
    The function will prompt for a session name after plotting the data within a VisIt window.
    
    """
#    dir of data
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    #Get end time variable
    i = 0
    t = []
    path = os.path.exists("batch/" + name + "_batch_%d.vts" % i)
    while path == True:
        t.append(i)
        i+=skip
        path = os.path.exists("batch/" + name + "_batch_%d.vts" % i)
    
    #Make dir for storing images
    if not os.path.exists("images"):
            os.makedirs("images")
    
    #Set the image dir
    img_dir = work_dir + "/images"
    
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
    visit.Launch(vdir=visual.visit_bin)
    
    print("Please orientate the Plot how you would like the image sequence to appear, then save the session file")
        
    # If a max min file exists assign a value if not set Max,Min = False
    if not os.path.exists("max_min_" + name + '.txt'):
        max,min = False, False
    else:
        # Import the max and min values for the data
        mm_array = np.loadtxt('max_min_' + name + '.txt')
        max = mm_array[0]
        min = mm_array[1]
    
    print('max',max)
    print('min',min)
        
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
        
    # View the VTK files so that the user can reorientate the Plot
    session_path,session_name = visual.view_vtk(work_dir,name,max,min)
    print(session_path)
       
    use_max = int(raw_input('Use max from max input file? (0 for No, 1 for Yes): '))
    use_min = int(raw_input('Use min from min input file? (0 for No, 1 for Yes): '))
    
    # Max min
    if use_max == 0:
        max = False
        min = False
    if use_min == 0:
        max = False
        min = False
    if (use_max == 1) or (use_min == 1):
        # Import the max and min values for the data
        mm_array = np.loadtxt('max_min_' + name +'.txt')
    if use_max == 1:
        max = mm_array[0]
    if use_min == 1:
        min = mm_array[1]  
    
    # Export an image sequence of the orientated data
    visual.draw_vtk(session_path,img_dir,name,t,session_name,max,min,skip)
    
    
#Vector draws a vector variable
def vector(name , skip = 1, path = None):
    
    """
    This fuction imports field parallel scalar data that has been converted into a vector, displays the data within a VisIt window and then renders out the image sequence.
    
    Inputs:
        name: name of variable converted
        skip: skip value used in conversion
        path: path to data if not located within working directory.
        
    Ouput:
        Image sequence of the VisIt plot.
    
    The function will prompt for a session name after plotting the data within a VisIt window.
    
    """
    #    dir of data
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    #Get end time variable
    i = 0
    t = []
    path = os.path.exists("batch/" + name + "_vector_%d.vts" % i)
    while path == True:
        t.append(i)
        i+=skip
        path = os.path.exists("batch/" + name + "_vector_%d.vts" % i)
    
    #Make dir for storing images
    if not os.path.exists("images"):
            os.makedirs("images")
    
    #Set the image dir
    img_dir = work_dir + "/images"
    
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
    visit.Launch(vdir=visual.visit_bin)
    
    print("Please orientate the Plot how you would like the image sequence to appear, then save the session file")
        
    # If a max min file exists assign a value if not set Max,Min = False
    if not os.path.exists("max_min_" + name + '.txt'):
        max,min = False, False
    else:
        # Import the max and min values for the data
        mm_array = np.loadtxt('max_min_' + name + '.txt')
        max = mm_array[0]
        min = mm_array[1]
    
    print('max',max)
    print('min',min)
        
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
        
    # View the VTK files so that the user can reorientate the Plot
    # Draw vtk file and let user orientate view and then save session file
    # returns the VisIt session name and location

    session_path,session_name = view_vector(work_dir,name,max,min)
    print(session_path)
       
    use_max = int(raw_input('Use max from max input file? (0 for No, 1 for Yes): '))
    use_min = int(raw_input('Use min from min input file? (0 for No, 1 for Yes): '))
    
    # Max min
    if use_max == 0:
        max = False
        min = False
    if use_min == 0:
        max = False
        min = False
    if (use_max == 1) or (use_min == 1):
        # Import the max and min values for the data
        mm_array = np.loadtxt('max_min_' + name +'.txt')
    if use_max == 1:
        max = mm_array[0]
    if use_min == 1:
        min = mm_array[1]  
    
    # Export an image sequence of the orientated data
    image_vector(session_path,img_dir,name,t,session_name,max,min,skip)
    
    
def eigen(name ,  path = None):
    """
    This fuction imports data from the eigensolver that has been converted in a coordinate system. Then the function displays the data within two VisIt windows, one for real data values and the other for imaginary data values, and then renders out the image sequence of both windows.
    
    The orientation of the real data is applied to the imaginary data, for the image sequence.
    
    Inputs:
        name: name of variable converted
        skip: skip value used in conversion
        path: path to data if not located within working directory.
        
    Ouput:
        Image sequence of the VisIt plot.
    
    The function will prompt for a session name after plotting the data within a VisIt window.
    
    
    """
    #    dir of data
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    
    #Get max time variable
    num_eig = len(visual.collect("t_array"))/2

    # Create list of the time values
    i = 0
    t = []
    while i < num_eig:
        path = os.path.exists("batch/" + name + "_eigen_%d.vts" % i)
        if path == True:
            t.append(i)
        i+=1
  
    #Make dir for storing images
    if not os.path.exists("images"):
            os.makedirs("images")
    
    #Set the image dir
    img_dir = work_dir + "/images"
    
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
    visit.Launch(vdir=visual.visit_bin)
    
    print("Please orientate the Plot how you would like the image sequence to appear, then save the session file")
    print("The orientation for the image sequence will be taken from window 1, so that both images sequences have the same orientation")

        
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
        
    # View the VTK files so that the user can reorientate the Plot
    session_path,session_name = visual.view_eigen(work_dir,name)
    print(session_path)
    
    # Export an image sequence of the orientated data
    visual.image_eigen(session_path,img_dir,name,t,session_name)
    
 
