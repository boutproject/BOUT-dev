import visual
import sys
import os
import visit
import numpy as np
import settings as set

# draw_vtk.py 

# Get the name of the variable to be plotted
#name = str(raw_input('Enter variable name: '))
name = str(sys.argv[1])

#Change working dir to Parent of VisBOUTIt
work_dir = os.getcwd()
os.chdir(work_dir + '/../')
work_dir = os.getcwd()

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

print set.visit_dir
#Launch VisIt
sys.path.insert(0,set.visit_dir)
import visit
visit.Launch(vdir=set.visit_bin)

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



# View the VTK files so that the user can reorientate the Plot
session_path,session_name = visual.view_vtk(work_dir,name,max,min)

print session_path
# Export an image sequence of the orientated data
visual.draw_vtk(session_path,img_dir,name,t,session_name,max,min)
