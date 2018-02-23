.. _sec-idl-routines:

IDL routines
============

List of IDL routines available in idllib. There are broadly three
categories of routine:

-  Completely general routines which could be useful outside BOUT++ work

   -  Data plotting and animation: **contour2** and **showdata**

   -  File reading and writing: **file\_open**, **file\_read** etc.

   -  User input and output: **get\_float**, **get\_integer**,
      **get\_yesno** and **str**

   -  FFT routines for integrating, differentiating and filtering:
      **fft\_integrate**, **fft\_deriv**, **fft\_filter**

-  Routines for BOUT++, but not specific to any application

   -  Modifying restart files: **expand\_restarts**, **scale\_restarts**
      and **split\_restarts**

   -  Processing 3D variables for input grid: **bout3dvar**

-  Routines specifically for tokamak simulations

   -  Reading A- and G-EQDSK format files into IDL: **read\_aeqdsk** and
      **read\_neqdsk**

   -  Plotting results: **polslice**, **plotpolslice**

Here the format is

**name**, arguments, [optional arguments]

-  | var = **bout3dvar** ( var )
   | Converts 3D variables to and from BOUT++’s Fourier representation
     which is used for input grids. By default converts from [x,y,z] to
     [x,y,f]

   -  **/reverse** Convert from [x,y,f] to [x,y,z]

   -  **nf**\ =nf Set number of frequencies in the result

   -  **nz**\ =nz When using /reverse, set number of Z points in the
      result

-  | var = **collect**\ ()
   | Read in data from a set of BOUT++ dump files

   -  **var** = “name of variable”

   -  **path** = “path/to/variable/”

   -  **xind**, **yind**, **zind**, **tind** = [min, max] index pairs

   -  **t\_array** = Output 1D array of times

-  | **contour2**, data [, x, y]
   | This is a replacement for the IDL contour which includes a scale
     color bar.

   -  **data** can be either 2D (x,y) or 3D (x,y,t). If data is 3D then
      the color is scaled to the entire range.

   -  **x** is an optional 2D (x,y) array of X coordinates

   -  **y** is an optional 2D (x,y) array of Y coordinates

   -  **t**\ =t is a time index for 3D data

   -  **nlev**\ =nlev

   -  **centre**\ =centre Make zero the middle of the color range (white
      if redblue)

   -  **redblue**\ =redblue Use a blue-white-red color scheme

   -  **revcolor**\ =revcolor Reverse color scheme

-  | **expand\_restarts**, newz
   | Increases the number of Z points in restart files. Together with
     scale\_restarts and split\_restarts, this makes it easier to modify
     a linear simulation as a start for non-linear runs.

   -  **newz** is the new value of NZ

   -  **path**\ =path Input path

   -  **output**\ =output Output path

   -  **format**\ =format File extension of output

-  | result = **fft\_deriv** ( var1d )
   | Calculates the derivative of a variable on a periodic domain.

-  result = **fft\_filter** (var, nf) Fourier filter a variable on a
   periodic domain. Arguments are a 1D variable and the number of
   Fourier components to keep

-  result = **fft\_integrate** ( var1d ) Integrates a variable on a
   periodic domain.

   -  **loop**\ =loop The loop integral is returned in this variable

-  | **file\_close**, handle
   | Close a file opened using file\_open()

-  | list = **file\_list** ( handle )
   | Return a list of variable names in the file

-  | integer = **file\_ndims** ( handle , “variable” )
   | Get the number of dimensions of a variable

-  | handle = **file\_open** ( “file” )
   | Open a NetCDF file.

   -  **/write** Open file for writing (default is read only)

   -  **/create** Create a new file, over-writing if already exists

-  var = **file\_read** ( handle, “variable” )

   -  **inds** = [xmin, xmax, ymin, ymax, ... ]

-  | float = **get\_float** ( “prompt” )
   | Ask the user for a float, using the given prompt

-  | integer = **get\_integer** ( “prompt” )
   | Ask the user for an integer

-  | integer = **get\_yesno** ( “prompt” )
   | Ask for a yes (1) or no (0) answer

-  | result = **gmres** ( x0, operator, b )
   | General Minimal Residual (GMRES)

   -  **x0** is the starting guess at the solution

   -  **operator**

   -  **b**

   Optional arguments

   -  **restart**\ =restart

   -  **max\_iter**\ =max\_iter

   -  **tol**\ =tol

   -  **stats**\ =stats

   -  **show**\ =show

   -  **output**\ =output

-  | result = **int\_func** ( [x,] f )
   | Integrate a function, always using the maximum number of
     grid-points possible for highest accuracy

-  | bool = **is\_pow2** ( value )
   | Returns 1 (true) if the given number is a power of 2, 0 (false)
     otherwise

-  | **plotpolslice**, var3d, grid
   | Takes a slice through a field-aligned tokamak domain, showing a
     poloidal cross-section.

   -  **var3d** is a 3D (x,y,z) variable to plot. Needs all of the
      points to work properly.

   -  **grid** is a structure from importing a grid file

   Optional arguments:

   -  **period**\ =period

   -  **zangle**\ =zangle

   -  **nlev**\ =nlev

   -  **yr**\ =yr

   -  **profile**\ =profile

   -  **output**\ =output

   -  **lines**\ =lines

   -  **linecol**\ =linecol

   -  **filter**\ =filter

-  | **polslice**, data, gridfile
   | Plots a 2D poloidal contour for single or double-null
     configurations, including color bar.

   -  **xstart**\ =xstart X index where the data begins. Useful if only
      part of the domain has been collected

   -  **ystart**\ =ystart Y index where data begins

-  | struct = **read\_aeqdsk** ( “filename” )
   | Reads an A-EQDSK file. Format is specified here:
     https://fusion.gat.com/THEORY/efit/a_eqdsk.html

-  | struct = **read\_neqdsk** ( “filename” )
   | Reads in an ’neqdsk’ or G-EQDSK formatted tokamak equilibrium file.
     Format of G-EQDSK file is specified here:
     https://fusion.gat.com/THEORY/efit/g_eqdsk.html

-  | stringarray = **regex\_extract** ( line, pattern )
   | Extract all matches to Regular Expression pattern contained in
     line. Useful for extracting numbers from FORTRAN-formatted text
     files.

   -  **line** Input string

   -  **pattern** Regular expression pattern to match

   -  **nmatch**\ =nmatch

-  | var = **reverse\_inds** ( var )
   | Reverse array indices e.g. ``arr[t,z,y,x] -> arr[x,y,z,t]``. Works
     on up to 5 dimensional variables

-  | **safe\_colors**
   | Sets the color table to useful values for plotting.

   -  **/first** Sets the first 10 colors to specific values, otherwise
      sets last 7

-  **scale\_restarts**, factor

   -  **path**\ =path Path to the restart files (default is current
      directory ’.’)

   -  **format**\ =format Specify what the file format is, otherwise
      goes on the file name

-  | **showdata**, data
   | Display animations of 1D,2D and 3D data. Defaults:

   -  2D data Animate a line plot

   -  3D data Animate a surface plot

   -  4D data Animate a poloidal cross-section (tokamaks only)

   Optional arguments:

   -  **/addsym** For 2D data (1D plots), add symbols to mark data
      points

   -  **az**\ =angle Rotate surface plots

   -  **/bw** Make contour plots grey scale

   -  **chars**\ =size character size

   -  **/contour** For 3D input, show color contour plot

   -  **delay**\ =time Time delay between plots (default 0.2 seconds)

   -  **/noscale** By default, all plots are on the same scale. This
      changes the scale for each plot’s range

   -  **profile**\ =array Background profile. Data is 3D: profile is 1D
      (X). Data is 4D -> profile is 2D (X,Y)

   -  **yr**\ =[min,max] Y range

-  | result = **sign** ( var )
   | This returns +1 if the variable is :math:`> 0`, -1 otherwise

-  **spectrum**

-  | **split\_restarts**, [nxpe], nype
   | split restart files between a different number of processors

   -  **nxpe** is an optional argument giving the number of processors
      in the X direction

   -  **nype** is the number of processors in the Y direction

   -  **path**\ =path Input path

   -  **output**\ =output Output path

   -  **format**\ =format File extension of output

-  | string = **str** ( value )
   | Convert a value to a string with whitespace trimmed. Arrays are
     converted to a comma-separated list in brackets.

-  | result = **zfamp** ( var4d )
   | Given a 4D variable [x,y,z,t], returns the Fourier amplitudes in
     [x,y,f,t]

-  | var = **zshift** ( var, shift )
   | Shifts a variable in the Z direction, useful for mapping between
     field-aligned and orthogonal coordinates.

   -  **period**\ =period How many domains fit in :math:`2\pi`. Default
      is 1 (full torus)

