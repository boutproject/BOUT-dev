# Analysis routines for Python
(Note that this document is out of date)

## Dependencies
Following packages must be installed for these
routines to work:
* **NumPy**           for arrays
* **netcdf4-python**  for NetCDF file reading

## Available routines
* **boutdata/**         BOUT++ data reading package
  * **collect**   Collect data from BOUT++ dump files


## Examples
Reading data from dump files:

```
from boutdata import *
ni = collect("Ni")
```
