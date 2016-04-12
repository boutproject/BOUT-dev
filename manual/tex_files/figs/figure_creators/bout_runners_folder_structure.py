#!/usr/bin/env python
"""
Manually creation of folder structures
To be used with python 2

Dependencies:
-----------
pygraphviz

mmag@fysik.dtu.dk
"""

__authors__ = 'Michael Loeiten'
__email__   = 'mmag@fysik.dtu.dk'
__version__ = '1.0'
__date__    = '21.01.2016'

import pygraphviz as pgv

# Initialize
tree = pgv.AGraph()

# Appendable lists
files     = []
dead_ends = []

# Default node attributes
tree.node_attr['shape']='box'
tree.node_attr['style']='bold'

# Adding nodes and edges
l0 = 'project'
l1 = ['data', 'source\nfiles', 'driver.py']
# Append the files
files.append(l1[1])
files.append(l1[2])
# Add the boxes to the mother node
for box in l1:
    tree.add_edge(l0,box)

l2 = ['solver1', 'solver2',\
      'BOUT.inp', 'run_log.txt']
# Append the files
files.append(l2[2])
files.append(l2[3])
# Add the boxes to the mother node
for box in l2:
    tree.add_edge('data', box)
tree.add_edge('solver2', 'solver2/...')
# Append the dead_end
de = l2[1] + '/...'
dead_ends.append(de)

l3 = ['method1', 'method2', 'solver1/...']
for box in l3:
    tree.add_edge('solver1', box)
tree.add_edge('method2', 'method2/...')
# Append the dead_end
de = l3[2]
dead_ends.append(de)
de = l3[1] + '/...'
dead_ends.append(de)

l4 = ['nout\ntimestep1', 'nout\ntimestep2', 'method1/...']
for box in l4:
    tree.add_edge('method1', box)
tree.add_edge('nout\ntimestep2', 'nout\ntimestep2/...')
# Append the dead_end
de = l4[2]
dead_ends.append(de)
de = l4[1] + '/...'
dead_ends.append(de)

l5 = ['mesh1', 'mesh2', 'nout\ntimestep1/...']
for box in l5:
    tree.add_edge('nout\ntimestep1', box)
tree.add_edge('mesh2', 'mesh2/...')
# Append the dead_end
de = l5[2]
dead_ends.append(de)
de = l5[1] + '/...'
dead_ends.append(de)

l6 = ['additional1', 'additional2', 'mesh1/...']
for box in l6:
    tree.add_edge('mesh1', box)
tree.add_edge('additional2', 'additional2/...')
# Append the dead_end
de = l6[2]
dead_ends.append(de)
de = l6[1] + '/...'
dead_ends.append(de)

l7 = ['grid_file1', 'grid_file2', 'additional1/...']
for box in l7:
    tree.add_edge('additional1', box)
tree.add_edge('grid_file2', 'grid_file2/...')
# Append the dead_end
de = l7[2]
dead_ends.append(de)
de = l7[1] + '/...'
dead_ends.append(de)

l8 = ['BOUT.inp\n(copy)', 'BOUT.log', 'BOUT.dmp',\
      'BOUT.restart', '(source_files\n(copy))', '(grid_file\n(copy))']
# Add l8 to the files list
for cur_file in l8:
    files.append(cur_file)
# Append them to the mother node
for box in l8:
    tree.add_edge('grid_file1', box)


# Change colors for the files
for the_file in files:
    member=tree.get_node(the_file)
#    member.attr['fontcolor'] = 'limegreen'
    member.attr['color'] = 'limegreen'

# Change colors for the dead_ends
for dead_end in dead_ends:
    member=tree.get_node(dead_end)
#    member.attr['fontcolor'] = 'darksalmon'
    member.attr['color'] = 'darksalmon'


# Print the graph
print(tree.string())

# Set layout
tree.layout('dot')
# Write to file
tree.draw('folder_tree.svg')
