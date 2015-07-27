# image.py
# Create image seqeuence from VisIt Session File
import visual
import sys
import os 
import numpy as np
#import ConfigParser as cp

def image():
    name = str(raw_input('Please enter variable name: '))
    session_name = str(raw_input('Please Enter session name: '))
    use_max = int(raw_input('Use max from max input file? (0 for No, 1 for Yes): '))
    use_min = int(raw_input('Use min from min input file? (0 for No, 1 for Yes): '))

    #Get the working directory
    work_dir = os.getcwd()
    #Get end time variable
    i = 0
    
    path = os.path.exists(work_dir + "/batch/" + name + "_batch_%d.vts" % i)
    while path == True:
        i +=1
        path = os.path.exists(work_dir + "/batch/" + name + "_batch_%d.vts" % i)
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

    session_path = work_dir + '/' + session_name + '.session'
    # Export an image sequence of the orientated data
    visual.draw_vtk(session_path,img_dir,name,t,session_name,max,min)

