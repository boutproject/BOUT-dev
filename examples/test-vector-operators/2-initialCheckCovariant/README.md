# 2-initialCheckCovariant

Contains an implementation which shows how to set up calculations using a
vector field, and how to plot these. Note that the routines here gives
identical result as in [1-initialCheckContravariant/](../1-initialCheckContravariant/)

Folders present are:

* [calculations](./calculations/) contains ipython notebooksof calculations of
  the gradient in cylindrical coordinates, and the mapping between
  contravariant cylindrical vectors and cartesian vectors.
    * The files can be viewed on the github page
    * If the equations does not render properly, try copying the url of the
      file to http://nbviewer.ipython.org/
* [cartesian](./cartesian/) contains the input file using cartesian geometry
* [cylindrical](./cylindrical/) contains the input file using cylindrical geometry
* [pythonRoutines](./pythonRoutines/) contains post-processing routines

Driver files present are:

* [driverCartesian.py](driverCartesian.py) bout_runners driver for the cartesian
  case
* [driverCylindrical.py](driverCylindrical.py) bout_runners driver for the
  cylindrical case
