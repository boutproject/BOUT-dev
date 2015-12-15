# 1-initialCheckContravariant

Contains an implementation which shows how to set up calculations using a
vector field, and how to plot these. Note that the routines here gives
identical result as in [2-initialCheckCovariant/](../2-initialCheckCovariant/)

Folders present are:

* [calculations](/calculations/) contains ipython notebooks (see
  [jupyter](http://jupyter.org/) for more details (can by the way be viewed on
  the github page)) of calculations of the gradient in cylindrical coordinates,
  and the mapping between contravariant cylindrical vectors and cartesian
  vectors.
* [cartesian](/cartesian/) contains the input file using cartesian geometry
* [cylindrical](/cylindrical/) contains the input file using cylindrical geometry
* [pythonRoutines](/pythonRoutines/) contains post-processing routines

Driver files present are:

* [driverCartesian.py](driverCartesian.py) bout_runners driver for the cartesian
  case
* [driverCylindrical.py](driverCylindrical.py) bout_runners driver for the
  cylindrical case
