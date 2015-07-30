# draw.py 
# Create image sequence from batch files
import visual
import sys
import os
import numpy as np

def draw(path = None):
#    dir of data
    if path == None:
        work_dir = os.getcwd()
    if path != None:
        work_dir = os.chdir(path)
        work_dir = os.getcwd()
    # Get the name of the variable to be plotted
    name = str(raw_input('Enter variable name: '))
    
    #Get end time variable
    i = 0
    path = os.path.exists("batch/" + name + "_batch_%d.vts" % i)
    while path == True:
        i +=1
        path = os.path.exists("batch/" + name + "_batch_%d.vts" % i)
    
    t = i
    #Make dir for storing images
    if not os.path.exists("images"):
            os.makedirs("images")
    
    #Set the image dir
    img_dir = work_dir + "/images"
    
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
    visit.Launch(vdir=visual.visit_bin)
    
    print "Please orientate the Plot how you would like the image sequence to appear, then save the session file"
        
    # If a max min file exists assign a value if not set Max,Min = False
    if not os.path.exists("max_min_" + name + '.txt'):
        max,min = False, False
    else:
        # Import the max and min values for the data
        mm_array = np.loadtxt('max_min_' + name + '.txt')
        max = mm_array[0]
        min = mm_array[1]
    
    print 'max',max
    print 'min',min    
    
    
    #Launch VisIt
    sys.path.insert(0,visual.visit_dir)
    import visit
    
    
    # View the VTK files so that the user can reorientate the Plot
    session_path,session_name = visual.view_vtk(work_dir,name,max,min)
    print session_path
       
    use_max = int(raw_input('Use max from max input file? (0 for No, 1 for Yes): '))
    use_min = int(raw_input('Use min from min input file? (0 for No, 1 for Yes): '))
    
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
              
    
    # Export an image sequence of the orientated data
    visual.draw_vtk(session_path,img_dir,name,t,session_name,max,min)
