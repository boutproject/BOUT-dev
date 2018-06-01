"""
image.py
Export an image seqeuence from VisIt Session File for scalar, vector and eigensolver data.
"""

 
from . import visual
import sys
import os 
import numpy as np

def image(name, skip = 1 , path = None):
    """
    This function imports a VisIt session of scalar data (with a single window) and then renders out the image sequence of the avaliable data
    
    Inputs:
        name: variable name of the converted data.
        path: path to the raw data, if not present within the current working directory.
        
    Output:
        An image sequence of the VisIt session.
        
    The Script prompts for session name and options for using the maximum and minimum values, for a Pseudocolor plot. The script then renders out the image sequence.
    
    """
    session_name = str(raw_input('Please Enter session name: '))
    use_max = int(raw_input('Use max from max input file? (0 for No, 1 for Yes): '))
    use_min = int(raw_input('Use min from min input file? (0 for No, 1 for Yes): '))

    #Get the working directory
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
        
    #Get end time variable
    i = 0
    t = []
    path = os.path.exists(work_dir + "/batch/" + name + "_batch_%d.vts" % i)
    while path == True:
        t.append(i)
        i +=skip
        path = os.path.exists(work_dir + "/batch/" + name + "_batch_%d.vts" % i)
    
    # Make dir for storing images
    if not os.path.exists("images"):
        os.makedirs("images")

    # Set the image dir
    img_dir = work_dir + "/images"
    
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
    visit.Launch(vdir=visual.visit_bin)
    
    #Max min
    if use_max == 0:
        max = False
        min = False
    if use_min == 0:
        max = False
        min = False
    if (use_max == 1) or (use_min == 1):
        #Import the max and min values for the data
        mm_array = np.loadtxt('max_min_' + name +'.txt')
    if use_max == 1:
        max = mm_array[0]
    if use_min == 1:
        min = mm_array[1]

    # Launch Visit Session
    session_path = work_dir + '/' + session_name + '.session'
    # Export an image sequence of data
    visual.draw_vtk(session_path,img_dir,name,t,session_name,max,min,skip)
    
def vector(name, skip = 1 , path = None):
    """
    This function imports a VisIt session of vector data (with a single window) and then renders out the image sequence of the avaliable data
    
    Inputs:
        name: variable name of the converted data.
        path: path to the raw data, if not present within the current working directory.
        
    Output:
        An image sequence of the VisIt session.
        
    The Script prompts for session name and options for using the maximum and minimum values, for a Pseudocolor plot. The script then renders out the image sequence.
    
    """
    session_name = str(raw_input('Please Enter session name: '))
    use_max = int(raw_input('Use max from max input file? (0 for No, 1 for Yes): '))
    use_min = int(raw_input('Use min from min input file? (0 for No, 1 for Yes): '))

    #Get the working directory
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
        
    #Get end time variable
    i = 0
    t = []
    path = os.path.exists(work_dir + "/batch/" + name + "_vector_%d.vts" % i)
    while path == True:
        t.append(i)
        i +=skip
        path = os.path.exists(work_dir + "/batch/" + name + "_vector_%d.vts" % i)
    
        #Max min
    if use_max == 0:
        max = False
        min = False
    if use_min == 0:
        max = False
        min = False
    if (use_max == 1) or (use_min == 1):
        #Import the max and min values for the data
        mm_array = np.loadtxt('max_min_' + name +'.txt')
    if use_max == 1:
        max = mm_array[0]
    if use_min == 1:
        min = mm_array[1]    
    
    # Make dir for storing images
    if not os.path.exists("images"):
        os.makedirs("images")

    # Create index of number of time slices (used for keeping time in the filenames and navigating in VisIt)
    indicies = np.arange(0,len(t))
        
    # Launch visit
    sys.path.insert(0,visual.visit_dir)
    import visit
    # Load session and initialise at time 0
    visit.Launch(vdir=visual.visit_bin)

    # Session path    
    session_path = work_dir + '/' + session_name + '.session'    
    

    # Set the image dir
    img_dir = work_dir + "/images" 

    # Make image session dir    
    outputdir = img_dir + '/' + session_name
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)    
    
    visit.RestoreSession(session_path,0)
    visit.SetTimeSliderState(0)
    # Export an Image sequence of the variable for every time base
    for i in indicies:
        time = i * skip
        visit.SetTimeSliderState(i) # Change timer slider
        VectorAtts = visit.VectorAttributes()
        # If user would like fixed max and mind then assign max and min values
        if max != 0:
            VectorAtts.max = max
            VectorAtts.maxFlag = 1
        if max == 0:
            VectorAtts.maxFlag = 0
        if min != 0:
            VectorAtts.min = min
            VectorAtts.minFlag = 1
        if min == 0:
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
    
# Create image sequence for eigen (2 window interface)
def eigen(name, path = None):
    """
    This function imports a VisIt session of eigensolver data with two windows (window 1: real data, window 2: imaginary data) and then renders out the image sequence of the avaliable data
    
    Inputs:
        name: variable name of the converted data.
        path: path to the raw data, if not present within the current working directory.
        
    Output:
        An image sequence of the VisIt session.
        
    The Script prompts for session name and options for using the maximum and minimum values, for a Pseudocolor plot. The script then renders out the image sequence.
    
    """
    session_name = str(raw_input('Please Enter session name: '))

    #Get the working directory
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
    
    # Make dir for storing images
    if not os.path.exists("images"):
        os.makedirs("images")

    # Set the image dir
    img_dir = work_dir + "/images"
    
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
    visit.Launch(vdir=visual.visit_bin)

    # Launch Visit Session
    session_path = work_dir + '/' + session_name + '.session'
    # Export an image sequence of data
    visual.image_eigen(session_path,img_dir,name,t,session_name)
