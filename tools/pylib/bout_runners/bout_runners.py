#!/usr/bin/env python

# NOTE: THE BOUT-RUNNERS ARE UNDER REVISION. THINGS MAY NOT WORK AS
#       EXPECTED

"""Classes for running several mpi-runs with BOUT++ at once.
   Read the docstring of 'basic_runner', or refer to the user manual of
   BOUT++ for more info. Examples can be found in
   BOUT/examples/bout_runners_example."""

# NOTE: This document uses folding. A hash-symbol followed by three {'s
# denotes the start of a fold, and a hash-symbol followed by three }'s
# denotes the end of a fold
__authors__ = 'Michael Loeiten'
__email__   = 'mmag@fysik.dtu.dk'
__version__ = '0.92beta'
__date__    = '01.07.2015'

# TODO: DELETE PLOTTING ROUTINES? BOUTUTILS?
# TODO: Check that nx, ny and nz not input when using grid file

# TODO: Check if you can delete these
#from __future__ import print_function
#from builtins import zip
#from builtins import str
#from builtins import range
#from builtins import object

import os
import re
import itertools
import glob
import timeit
import datetime
from numbers import Number
import numpy as np
from subprocess import check_output
from boututils import shell, launch, getmpirun
from boututils.datafile import DataFile

# TODO: Things to add in __init__ and appropriate places
#       The following can already be set with the "additional" member
#       data
#       1.  Add all time integration solver option as given in "Table 6: Time
#           integration solver options" in the user manual
#       2.  Add laplace solver options

# TODO: Add additional error checks
#       1. For the branch cuts (to check if correct NXPE and NYPE)

# TODO: Make a check to see if the following is in the grid file
#       1. Check if branch cut is in the grid file

# TODO: Move allow_size_modification = False to execute_run, as it is
#       better to also write the run_no there?

# TODO: Test for python 2

# TODO: Rewrite documentation in the manual (also the nice chart)

# TODO: Make additional without combo

# FIXME: qsub does not always delete the clean-up files??
#        Fixed for basic qsub (test it), fix for the rest
# TODO: Make it possible to give a function to the waiting routine in
#       qsub runners, so that it is possible to give a function which
#       will be run when a job has completed
# TODO: Make qsub usable on different clusters (and update documentation)
#       Can be done by checking the current cluster? (Need to set the
#       path to the correct libraries, and use the correct MPI runner)
# TODO: Submit and postprocess to a different queue
# TODO: Check if it is possible to use qsub with dependencies (only
#       found depricated documentation)

#{{{class basic_runner
# As a child class uses the super function, the class must allow an
# object as input
class basic_runner(object):
#{{{docstring
    """Class for mpi running one or several runs with BOUT++.
    Calling self.execute_runs() will run your BOUT++ program with the possible
    combinations given in the member data using the mpi runner.

    Before each run, a folder system, based on the member data, rooted
    in self.__directory, will be created. The BOUT.inp of self.__directory
    is then copied to the execution folder.

    A log-file for the run is stored in self.__directory

    By default self.__directory = 'data', self.__nproc = 1 and
    self.__allow_size_modification = False

    self.__program_name is by default set to the same name as any .o files in the
    folder where an instance of the object is created. If none is found
    the creator tries to run make.

    All other data members are set to None by default.

    The data members will override the corresponding options given in
    self.__directory/BOUT.inp.

    See the doctring of the constructor (__int__) for options.
    See BOUT/examples/bout_runners_example for examples."""
#}}}

#{{{__init__
    def __init__(self,\
                 nproc      = 1,\
                 directory  = 'data',\
                 solver     = None,\
                 mms        = None,\
                 atol       = None,\
                 rtol       = None,\
                 mxstep     = None,\
                 grid_file  = None,\
                 nx         = None,\
                 ny         = None,\
                 nz         = None,\
                 zperiod    = None,\
                 zmin       = None,\
                 zmax       = None,\
                 dx         = None,\
                 dy         = None,\
                 dz         = None,\
                 MXG        = None,\
                 MYG        = None,\
                 ixseps1    = None,\
                 ixseps2    = None,\
                 jyseps1_1  = None,\
                 jyseps1_2  = None,\
                 jyseps2_1  = None,\
                 jyseps2_2  = None,\
                 symGlobX   = None,\
                 symGlobY   = None,\
                 ddx_first  = None,\
                 ddx_second = None,\
                 ddx_upwind = None,\
                 ddx_flux   = None,\
                 ddy_first  = None,\
                 ddy_second = None,\
                 ddy_upwind = None,\
                 ddy_flux   = None,\
                 ddz_first  = None,\
                 ddz_second = None,\
                 ddz_upwind = None,\
                 ddz_flux   = None,\
                 nout       = None,\
                 timestep   = None,\
                 additional = None,\
                 series_add = None,\
                 restart    = None,\
                 cpy_source = None,\
                 cpy_grid   = None,\
                 sort_by    = None,\
                 make       = None,\
                 allow_size_modification = False):
        #{{{docstring
        """The constructor of the basic_runner.

        All the member data is set to None by default. If the
        datamembers are not set, the values from BOUT.inp will be used.
        The exeption is nproc (default = 1), directory (default =
        'data') and allow_size_modification (default = False), which
        always needs to be set.

        Input:
        nproc       -    The number of processors to use in the mpirun (int)
        directory   -    The directory of the BOUT.inp file (str)
        solver      -    The solver to be used in the runs (str or
                         iterable)
        mms         -    Whether or not mms should be run (bool)
        atol        -    Absolute tolerance (number or iterable)
        rtol        -    Relative tolerance (number or iterable)
        mxstep      -    Max internal step pr output step (int or
                         iterable)
        grid_file   -    The grid file (str or iterable)
        nx          -    Number of nx in the run (int or iterable)
        ny          -    Number of ny in the run (int or iterable)
        nz          -    Number of nz in the run (int or iterable)
        zperiod     -    Domain size in  multiple of fractions of 2*pi
                         (int or iterable)
        zmin        -    Minimum range of the z domain
        zmax        -    Maximum range of the z domain
        dx          -    Grid size in the x direction (Number or iterable)
        dy          -    Grid size in the x direction (Number or iterable)
        dz          -    Grid size in the x direction (Number or iterable)
        MXG         -    The number of guard cells in the x direction
                         (int)
        MYG         -    The number of guard cells in the y direction
                         (int)
        ixseps1     -    Separatrix location for 'upper' divertor (int
                         or iterable)
        ixseps2     -    Separatrix location for 'lower' divertor (int
                         or iterable)
        jyseps1_1   -    Branch cut location 1_1 [see user's manual for
                         details] (int or iterable)
        jyseps1_2   -    Branch cut location 1_2 [see user's manual for
                         details] (int or iterable)
        jyseps2_1   -    Branch cut location 2_1 [see user's manual for
                         details] (int or iterable)
        jyseps2_2   -    Branch cut location 2_2 [see user's manual for
                         details] (int or iterable)
        symGlobX    -    symmetricGLobalX: x defined symmetrically between
                         0 and 1 (bool)
        symGlobY    -    symmetricGLobalX: y defined symmetrically (bool)
        ddx_first   -    Method used for for first ddx terms (str or iterable)
        ddx_second  -    Method used for for second ddx terms (str or iterable)
        ddx_upwind  -    Method used for for upwind ddx terms (str or iterable)
        ddx_flux    -    Method used for for flux ddx terms (str or iterable)
        ddy_first   -    Method used for for first ddy terms (str or iterable)
        ddy_second  -    Method used for for second ddy terms (str or iterable)
        ddy_upwind  -    Method used for for upwind ddy terms (str or iterable)
        ddy_flux    -    Method used for for flux ddy terms (str or iterable)
        ddz_first   -    Method used for for first ddz terms (str or iterable)
        ddz_second  -    Method used for for second ddz terms (str or iterable)
        ddz_upwind  -    Method used for for upwind ddz terms (str or iterable)
        ddz_flux    -    Method used for for flux ddz terms (str or iterable)
        nout        -    Number of outputs stored in the *.dmp.* files
                         (int or iterable)
        timestep    -    The time between each output stored in the
                         *.dmp.* files (int or iterable)
        additional  -    Additional option for the run given on the form
                         ('section_name','variable name', values) or as
                         iterable on the same form, where values can be
                         any value or string or an iterable of those
        series_add  -    The same as above, with the exception that
                         no combination will be performed between the elements
                         during a run
        restart     -    Wheter or not to use the restart files
                         ('overwrite' or 'append')
        cpy_source  -    Wheter or not to copy the source files to the
                         folder of the *.dmp.* files (bool)
        cpy_grid    -    Wheter or not to copy the grid files to the
                         folder of the *.dmp.* files (bool)
        sort_by     -    Defining what will be the fastest running
                         variable in the run, which can be useful if one
                         for example would like to 'group' the runs before
                         sending it to a post processing function (see
                         the docstring of the run function for more
                         info). The possibilities are
                        'spatial_domain',
                        'temporal_domain',
                        'solver',
                        'ddx_first',
                        'ddx_second',
                        'ddx_upwind',
                        'ddx_flux',
                        'ddy_first',
                        'ddy_second',
                        'ddy_upwind',
                        'ddy_flux',
                        'ddz_first',
                        'ddz_second',
                        'ddz_upwind',
                        'ddz_flux',
                        any 'variable_name' from additional or series_add
                        an iterable consisting of several of these. If
                        an iterable is given, then the first element is
                        going to be the fastest varying variable, the
                        second element is going to be the second fastest
                        varying variable and so on.
        make        -   Whether or not to make the progam (bool)

        allow_size_modification - Whether or not to allow bout_runners
                                  modify nx and ny in order to find a
                                  valid split of the domain (bool)

        Please note:
        - The length of nx, ny and nz must be the same
        - The length of timestep and nout must be the same
        - If nx, ny or nz is given in the grid file, bout_runners will use
          these values.
        """
        #}}}

        # Setting the member data
        self.__nproc      = nproc
        self.__directory  = directory
        self.__solver     = self.__set_member_data(solver)
        self.__mms        = mms
        self.__atol       = self.__set_member_data(atol)
        self.__rtol       = self.__set_member_data(rtol)
        self.__mxstep     = self.__set_member_data(mxstep)
        self.__grid_file  = self.__set_member_data(grid_file)
        self.__nx         = self.__set_member_data(nx)
        self.__ny         = self.__set_member_data(ny)
        self.__nz         = self.__set_member_data(nz)
        self.__zperiod    = self.__set_member_data(zperiod)
        self.__zmin       = self.__set_member_data(zmin)
        self.__zmax       = self.__set_member_data(zmax)
        self.__dx         = self.__set_member_data(dx)
        self.__dy         = self.__set_member_data(dy)
        self.__dz         = self.__set_member_data(dz)
        self.__MXG        = MXG
        self.__MYG        = MYG
        self.__ixseps1    = self.__set_member_data(ixseps1)
        self.__ixseps2    = self.__set_member_data(ixseps2)
        self.__jyseps1_1  = self.__set_member_data(jyseps1_1)
        self.__jyseps1_2  = self.__set_member_data(jyseps1_2)
        self.__jyseps2_1  = self.__set_member_data(jyseps2_1)
        self.__jyseps2_2  = self.__set_member_data(jyseps2_2)
        self.__symGlobX   = symGlobX
        self.__symGlobY   = symGlobY
        self.__ddx_first  = self.__set_member_data(ddx_first)
        self.__ddx_second = self.__set_member_data(ddx_second)
        self.__ddx_upwind = self.__set_member_data(ddx_upwind)
        self.__ddx_flux   = self.__set_member_data(ddx_flux)
        self.__ddy_first  = self.__set_member_data(ddy_first)
        self.__ddy_second = self.__set_member_data(ddy_second)
        self.__ddy_upwind = self.__set_member_data(ddy_upwind)
        self.__ddy_flux   = self.__set_member_data(ddy_flux)
        self.__ddz_first  = self.__set_member_data(ddz_first)
        self.__ddz_second = self.__set_member_data(ddz_second)
        self.__ddz_upwind = self.__set_member_data(ddz_upwind)
        self.__ddz_flux   = self.__set_member_data(ddz_flux)
        self.__nout       = self.__set_member_data(nout)
        self.__timestep   = self.__set_member_data(timestep)
        self.__additional = additional
        self.__series_add = series_add
        self.__restart    = restart
        self.__cpy_source = cpy_source
        self.__cpy_grid   = cpy_grid
        self.__sort_by    = self.__set_member_data(sort_by)
        self.__make       = make
        self.__allow_size_modification = allow_size_modification

        # Make some space to distinguish from the rest of the terminal
        print("\n")

        # Initializing self.__warnings and self.__error
        # self.__warnings will be filled with warnings
        # self.__errors will be filled with errors
        # The warnings and errors will be printed when the destructor is called
        self.__warnings = []
        self.__errors   = []

        # Check if make is a boolean
        if self.__make != None:
            if type(self.__make) != bool:
                self.__errors.append("TypeError")
                raise TypeError("make must be boolean if set")

        # self.__additional must be on a special form (see
        # basic_error_checker).
        if self.__additional != None:
            if not(hasattr(self.__additional, "__iter__")) or\
               (type(self.__additional) == str) or\
               (type(self.__additional) == dict):
                # Put additional as a double iterable
                self.__additional = [(self.__additional)]
            else:
                if not(hasattr(self.__additional[0], "__iter__")) or\
                   (type(self.__additional[0]) == str) or\
                   (type(self.__additional) == dict):
                    # Put self.__additional as an iterable
                    self.__additional = [self.__additional]
        # Do the same for series_add
        if self.__series_add != None:
            if not(hasattr(self.__series_add, "__iter__")) or\
               (type(self.__series_add) == str) or\
               (type(self.__series_add) == dict):
                # Put series_add as a double iterable
                self.__series_add = [(self.__series_add)]
            else:
                if not(hasattr(self.__series_add[0], "__iter__")) or\
                   (type(self.__series_add[0]) == str) or\
                   (type(self.__series_add) == dict):
                    # Put self.__series_add as an iterable
                    self.__series_add = [self.__series_add]

        # Set self.__program_name from the *.o file. Make the program if
        # the *.o file is not found
        self.__set_program_name()

        # Make the file if make is True
        if self.__make:
            self.__run_make()

        # Obtain the MPIRUN
        self.__MPIRUN = getmpirun()

        # The run type is going to be written in the run.log file
        self.__run_type = 'basic'

        # Check if the instance is set correctly
        self.__check_for_basic_instance_error()
#}}}

#{{{__del__
    def __del__(self):
        """The destructor will print all the warning and error messages"""

        # Switch to see if error occured
        error_occured = False

        # If errors occured
        if len(self.__errors) > 0:
            message = "! A " + self.__errors[0] + " occured. !"
            # Find the boarder length
            len_boarder = len(message)
            # Print the message
            print("\n"*2 + "!"*len_boarder)
            print(message)
            print('!'*len_boarder + "\n"*2)
            error_occured = True
        if len(self.__warnings) > 0:
            print('\n'*3 + 'The following WARNINGS were detected:')
            print('-'*80)
            for warning in self.__warnings:
                print(warning + '\n')
            print('-'*80 + '\n'*3)
        elif len(self.__warnings) > 0 and not(error_occured):
            print('\n'*3 + ' ' + '~'*69)
            print("| No WARNINGS detected before instance destruction in"+\
                  " 'bout_runners'. |")
#}}}

#{{{execute_run
    def execute_runs(self,\
                     remove_old = False,\
                     post_processing_function = None,\
                     post_process_after_every_run = True,\
                     **kwargs):
        #{{{docstring
        """
        Makes a run for each of the combination given by the member data.

        Input
        remove_old                  - boolean telling whether old run
                                      files should be deleted or not
        post_processing_function    - a function to be called after
                                      one or several run. This function
                                      must accept the string of
                                      self.__dmp_folder if
                                      post_process_after_each_run is
                                      True, and a list of dmp folders if
                                      post_process_after_each_run is
                                      False
        post_process_after_each_run - boolean telling whether
                                      post_processing_function
                                      should be called after each run
                                      (if True), or after the number of
                                      runs decided by self.__sort_by
                                      (see the constructor of
                                      basic_runner for more info)
        **kwargs                    - parameters to be passed to the
                                      post_processing_function
        """
        #}}}

        # Check for errors in the run function
        self.__error_check_for_run_input(remove_old,\
                                         post_processing_function,\
                                         post_process_after_every_run)

        # Create the run log
        self.__create_run_log()

        # We check that the given combination of nx and ny is
        # possible to perform with the given nproc
        if (self.__nx != None) and (self.__ny != None):
            self.__get_correct_domain_split()

        # Get the combinations of the member functions
        possibilities = self.__get_possibilities()
        combinations = self.__get_combinations(possibilities)

        # If we are not running the post processing function after every
        # run, make an appendable list over all the runs which will be
        # passed as an input parameter to the post processing function
        if not(post_process_after_every_run):
            list_of_dmp_folders = []

        # Print either 'now running' or 'now submitting'
        self.__print_run_or_submit()

        # Set self.__len_group if post_processing_function is set, but
        # self.__sort_by is None
        if (post_processing_function != None) and\
           (not(post_process_after_every_run)) and\
           (self.__len_group == None):
               # self.__len_group is to a number by __get_swapped_input_list
               # (which is called if self.__sort_by != None)
               # If there are no sorting, self.__len_group will be None
               # We will make self.__len_group the length of the
               # number of runs here
               self.__len_group = len(combinations)

        # The run
        for run_no, combination in enumerate(combinations):

            # Get the folder to store the data
            self.__prepare_dmp_folder(combination)

            if remove_old:
                # Remove old data
               self.__remove_data()

            # Copy the grid (if any) if cpy_grid files is True We must
            # do this after self.__remove_data as this removes all *.nc
            # files
            if (self.__cpy_grid) and (self.__grid_file != None):
                combination_list = combination.split()
                # Loop through all the combinations
                for elem in combination_list:
                    # Find the grid
                    if 'grid' in elem:
                        # Remove grid=, so that only the path remains
                        cur_grid = elem.replace('grid=', '')
                        # Copy the grid file
                        command = 'cp ' + cur_grid + ' ' + self.__dmp_folder + '/'
                        shell(command)

            # Check if the run has been performed previously
            do_run = self.__check_if_run_already_performed()
            # Do the actual runs
            if do_run:
                # Call the driver for a run
                self.__run_driver(combination, run_no)

            # If we would like to call a post_processing function
            if post_processing_function != None:
                if post_process_after_every_run:
                    # Call the post processing function
                    self.__call_post_processing_function(\
                            function = post_processing_function,\
                            folders  = self.__dmp_folder,\
                            **kwargs)
                else:
                    # Append the dmp folder to the list of dmp folders
                    list_of_dmp_folders.append(self.__dmp_folder)
                    # If the run_no+1 is divisible by self.__len_group
                    if ((run_no+1) % self.__len_group == 0):
                        # Call the post processing function
                        self.__call_post_processing_function(\
                                function = post_processing_function,\
                                folders  = list_of_dmp_folders,\
                                **kwargs)
                        # Reset the list_of_dmp_folders
                        list_of_dmp_folders = []

        # post_run defines what to be done after the runs have finished/
        # been submitted (if anything)
        self.__post_run(**kwargs)
#}}}

#{{{__run_driver
    def __run_driver(self, combination, run_no):
        """The machinery which actually performs the run"""

        # Get the time when the run starts
        start = datetime.datetime.now()
        # Do the run
        output, run_time = self.__single_run(combination)
        # Print info to the log file for the runs
        self.__append_run_log(start, run_no, run_time)
        print('\n')
#}}}

#{{{ Functions called by the constructor
#{{{__set_member_data
    def __set_member_data(self, input_parameter):
        """Returns the input_parameter as a list if it is different than None,
        and if it is not iterable"""

       # If the input_data is not set, the value in BOUT.inp will
       # be used
        if input_parameter != None:
            # If the input_data is not an iterable, or if it is a
            # string: Put it to a list
            if not(hasattr(input_parameter, "__iter__")) or\
               (type(input_parameter)) == str:
                input_parameter = [input_parameter]

        return input_parameter
#}}}

#{{{__set_program_name
    def __set_program_name(self):
        """Set self.__program_name from the *.o file. Make the program if
        the *.o file is not found"""

        # Find the *.o file
        o_files = glob.glob("*.o")
        if len(o_files) > 0:
            # Pick the first instance as the name
            self.__program_name = o_files[0].replace('.o', '')
        else:
            # Check if there exists a make
            make_file = glob.glob("*make*")
            if len(make_file) > 0:
                # Run make
                self.__run_make()
                # Set the make flag to False, so it is not made again
                self.__make = False
                # Search for the .o file again
                o_files = glob.glob("*.o")
                if len(o_files) > 0:
                    self.__program_name = o_files[0].replace('.o', '')
                else:
                    self.__program_name = False
                    message = 'The constructor could not make your'+\
                              ' program'
                    self.__errors.append("RuntimeError")
                    raise RuntimeError(message)
            else:
                self.__errors.append("RuntimeError")
                raise RuntimeError("No make file found in current" +\
                                   " directory")
#}}}

#{{{__check_for_basic_instance_error
    def __check_for_basic_instance_error(self):
        """Check if there are any type errors when creating the object"""

        #{{{Check if nproc and directory has the correct type
        # nproc and directory is set by default, however, we must check that
        # the user has not given them as wrong input
        if type(self.__nproc) != int:
            message  = "nproc is of wrong type\n"+\
                       "nproc must be given as an int"
            self.__errors.append("TypeError")
            raise TypeError(message)
        if type(self.__directory) != str:
            message  = "directory is of wrong type\n"+\
                       "directory must be given as a str"
            self.__errors.append("TypeError")
            raise TypeError(message)
        #}}}

        #{{{Check if MXG and MYG has the correct type
        # Check if MXG and MYG is given as a single int
        # One should not be allowed to specify MXG and MYG as an
        # iterable, as MXG is used to find the correct split, and
        # because it in principle could be incompatible with the method
        # (first, second etc.) used
        check_if_int = [\
                        (self.__MXG, 'MXG'),\
                        (self.__MYG, 'MYG'),\
                       ]
        self.__check_for_correct_type(var = check_if_int,\
                                      the_type = int,\
                                      allow_iterable = False)
        #}}}

        #{{{Check if BOUT.inp exsists in the self.__directory
        # Check if there are any BOUT.inp files in the self.__directory
        inp_file = glob.glob(self.__directory + "/BOUT.inp")
        if len(inp_file) == 0:
            self.__errors.append("RuntimeError")
            raise RuntimeError("No BOUT.inp files found in '" +\
                                self.__directory + "'")
        #}}}

        #{{{Check if grid_file are strings, and that they exsist
        if self.__grid_file != None:
            # Set a variable which is has length over one if the test fails
            not_found = []
            if type(self.__grid_file) == str:
                # See if the grid_file can be found
                grid_file = glob.glob(self.__grid_file)
                # The grid_file cannot be found
                if len(grid_file) == 0:
                    not_found.append(self.__grid_file)
            # If several grid files are given
            elif hasattr(self.__grid_file, "__iter__"):
                for elem in self.__grid_file:
                    # See if the grid_file can be found
                    grid_file = glob.glob(elem)
                    # The grid_file cannot be found
                    if len(grid_file) == 0:
                        not_found.append(elem)
            if len(not_found) > 0:
                message =  "The following grid files were not found\n"
                message += "\n".join(not_found)
                self.__errors.append("RuntimeError")
                raise RuntimeError(message)
        #}}}

        #{{{Check nx, ny, nz, zperiod, nout, mxstep, separatrix are int/iterable
        check_if_int = [\
            (self.__nx       , 'nx')       ,\
            (self.__ny       , 'ny')       ,\
            (self.__nz       , 'nz')       ,\
            (self.__zperiod  , 'zperiod')  ,\
            (self.__nout     , 'nout')     ,\
            (self.__mxstep   , 'mxstep')   ,\
            (self.__ixseps1  , 'ixseps1')  ,\
            (self.__ixseps2  , 'ixseps2')  ,\
            (self.__jyseps1_1, 'jyseps1_1'),\
            (self.__jyseps1_2, 'jyseps1_2'),\
            (self.__jyseps2_1, 'jyseps2_1'),\
            (self.__jyseps2_2, 'jyseps2_2'),\
            ]

        self.__check_for_correct_type(var = check_if_int,\
                                      the_type = int,\
                                      allow_iterable = True)
        #}}}

        #{{{Check timestep, atol, rtol, zmin/max, dx, dy, dz is Number/iterable
        # Check if the following is a number
        check_if_number = [\
            (self.__timestep, 'timestep'),\
            (self.__zmin    , 'zmin')    ,\
            (self.__zmax    , 'zmax')    ,\
            (self.__dx      , 'dx')      ,\
            (self.__dy      , 'dy')      ,\
            (self.__dz      , 'dz')      ,\
            (self.__atol    , 'atol')    ,\
            (self.__rtol    , 'rtol')     \
            ]

        self.__check_for_correct_type(var = check_if_number,\
                                    the_type = Number,\
                                    allow_iterable = True)
        #}}}

        #{{{Check if solver, grid_file, methods and sort_by is str/list of str
        # Check if instance is string, or an iterable containing strings
        check_if_string = [\
            (self.__solver    , 'solver')    ,\
            (self.__grid_file , 'grid_file') ,\
            (self.__ddx_first , 'ddx_first') ,\
            (self.__ddx_second, 'ddx_second'),\
            (self.__ddx_upwind, 'ddx_upwind'),\
            (self.__ddx_flux  , 'ddx_flux')  ,\
            (self.__ddy_first , 'ddy_first') ,\
            (self.__ddy_second, 'ddy_second'),\
            (self.__ddy_upwind, 'ddy_upwind'),\
            (self.__ddy_flux  , 'ddy_flux')  ,\
            (self.__ddz_first , 'ddz_first') ,\
            (self.__ddz_second, 'ddz_second'),\
            (self.__ddz_upwind, 'ddz_upwind'),\
            (self.__ddz_flux  , 'ddz_flux')  ,\
            (self.__sort_by   , 'sort_by')    \
            ]

        self.__check_for_correct_type(var = check_if_string,\
                                    the_type = str,\
                                    allow_iterable = True)
        #}}}

        #{{{Check if solver is set to the correct possibility
        # Check if the solver is possible
        # From /include/bout/solver.hxx
        possible_solvers = [\
            'cvode',\
            'pvode',\
            'ida',\
            'petsc',\
            'karniadakis',\
            'rk4',\
            'euler',\
            'rk3ssp',\
            'power',\
            'arkode'\
            ]

        # Do the check if the solver is set
        if self.__solver != None:
            self.__check_if_set_correctly(var = (self.__solver, 'solver'),\
                                          possibilities = possible_solvers)
        #}}}

        #{{{Check if the methods is set to the correct possibility
        # Check if ddx or ddy is possible
        possible_method = [\
            'C2',\
            'C4',\
            'W2',\
            'W3'\
            ]

        # Make a list of the variables
        the_vars = [\
            (self.__ddx_first , 'ddx_first') ,\
            (self.__ddx_second, 'ddx_second'),\
            (self.__ddy_first , 'ddy_first') ,\
            (self.__ddy_second, 'ddy_second')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] != None:
                self.__check_if_set_correctly(var           = var,\
                                              possibilities = possible_method)

        # Check if ddz is possible
        possible_method.append('FFT')

        # Make a list of the variables
        the_vars = [\
            (self.__ddz_first , 'ddz_first') ,\
            (self.__ddz_second, 'ddz_second') \
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] != None:
                self.__check_if_set_correctly(var           = var,\
                                              possibilities = possible_method)

        # Check for upwind terms
        possible_method = [\
            'U1',\
            'U4',\
            'W3'\
            ]

        # Make a list of the variables
        the_vars = [\
            (self.__ddx_upwind, 'ddx_upwind'),\
            (self.__ddy_upwind, 'ddy_upwind'),\
            (self.__ddz_upwind, 'ddz_upwind')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] != None:
                self.__check_if_set_correctly(var          = var,\
                                              possibilities = possible_method)

        # Check for flux terms
        possible_method = [\
            'SPLIT',\
            'NND'\
            ]

        # Make a list of the variables
        the_vars = [\
            (self.__ddx_flux  , 'ddx_flux'),\
            (self.__ddy_flux  , 'ddy_flux'),\
            (self.__ddz_flux  , 'ddz_flux')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] != None:
                self.__check_if_set_correctly(var           = var,\
                                              possibilities = possible_method)
        #}}}

        #{{{Check if sort_by is set to the correct possibility
        possible_sort_by = [\
            'spatial_domain',\
            'temporal_domain',\
            'solver',\
            'ddx_first',\
            'ddx_second',\
            'ddx_upwind',\
            'ddx_flux',\
            'ddy_first',\
            'ddy_second',\
            'ddy_upwind',\
            'ddy_flux',\
            'ddz_first',\
            'ddz_second',\
            'ddz_upwind',\
            'ddz_flux',\
            ]

        # Append the additionals and series_add
        # If additional is set
        if self.__additional != None:
            for additional in self.__additional:
                # The additional now contains a tuple of three elements
                # We would like to extract the name of them, and append
                # it to the possibilities list
                possible_sort_by.append(additional[0])
        # Do the same for series_add
        if self.__series_add != None:
            for series_add in self.__series_add:
                # The series_add now contains a tuple of three elements
                # We would like to extract the name of them, and append
                # it to the possibilities list
                possible_sort_by.append(series_add[0])

        # Make a list of the variables
        the_vars = [\
            (self.__sort_by, 'sort_by')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] != None:
                self.__check_if_set_correctly(var           = var,\
                                              possibilities = possible_sort_by)
        #}}}

        #{{{Check if restart is set correctly
        if self.__restart != None:
            if type(self.__restart) != str:
                self.__errors.append("TypeError")
                raise TypeError ("restart must be set as a string when set")

        possible_method = [\
            'overwrite',\
            'append'\
            ]

        # Make a list of the variables
        the_vars = [\
            (self.__restart, 'restart')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] != None:
                self.__check_if_set_correctly(var           = var,\
                                              possibilities = possible_method)
        #}}}

        #{{{Check for options set in both member data and in the grid file
        if self.__grid_file != None:
            # Check if the following variables are found in the grid
            # file
            check_if_in_grid =[\
                    (self.__nx       , "nx")               ,\
                    (self.__ny       , "ny")               ,\
                    (self.__nz       , "nz")               ,\
                    (self.__dx       , "dx")               ,\
                    (self.__dy       , "dy")               ,\
                    (self.__dz       , "dz")               ,\
                    (self.__MXG      , "MXG")              ,\
                    (self.__MYG      , "MYG")              ,\
                    (self.__ixseps1  , "ixseps1")          ,\
                    (self.__ixseps2  , "ixseps2")          ,\
                    (self.__jyseps1_1, "jyseps1_1")        ,\
                    (self.__jyseps1_2, "jyseps1_2")        ,\
                    (self.__jyseps2_1, "jyseps2_1")        ,\
                    (self.__jyseps2_2, "jyseps2_2")        ,\
                    (self.__symGlobX , "symmmetricGlobalX"),\
                    (self.__symGlobY , "symmmetricGlobalY") \
                    ]
            for var in check_if_in_grid:
                # If the variable is set
                if var[0] != None:
                    # Loop through the grid files
                    for grid_file in self.__grid_file:
                        # Open (and automatically close) the grid files
                        f = DataFile(grid_file)
                        # Search for mesh data in the grid file
                        grid_variable = f.read(var[1])
                        # If the variable is found
                        if grid_variable != None:
                            self.__errors.append("TypeError")
                            message  = var[1] + " was specified both in the "
                            message += "driver and in the grid file.\n"
                            message += "Please remove " + var[1]
                            message += " from the driver if you would "
                            message += "like to run with a grid file."
                            raise TypeError(message)
        #}}}

        #{{{If grid files are set: Use nx, ny and nz values in the grid file
        if self.__grid_file != None:
            # Make a dict of appendable lists
            spatial_domain = {'nx':[], 'ny':[], 'nz':[]}
            for grid_file in self.__grid_file:
                # Open (and automatically close) the grid files
                f = DataFile(grid_file)
                # Search for nx, ny and nz in the grid file
                mesh_types = ["nx", "ny", "nz"]
                for mesh_type in mesh_types:
                    grid_variable = f.read(mesh_type)
                    # If the variable is found
                    if grid_variable != None:
                        spatial_domain[mesh_type].append(grid_variable)
            # Check that the lengths of nx, ny and nz are the same
            # unless they are not found
            len_nx = len(spatial_domain['nx'])
            len_ny = len(spatial_domain['ny'])
            len_nz = len(spatial_domain['nz'])
            if len_nx != 0:
                self.__nx = spatial_domain['nx']
            if len_ny != 0:
                self.__ny = spatial_domain['ny']
            if len_nz != 0:
                self.__nz = spatial_domain['nz']
        #}}}

        #{{{Check that nx, ny and nz are of the same length
        if self.__nx != None and self.__ny != None:
            self.__check_if_same_len((self.__nx, 'nx'), (self.__ny, 'ny'))
        if self.__nx != None and self.__nz != None:
            self.__check_if_same_len((self.__nx, 'nx'), (self.__nz, 'nz'))
        if self.__ny != None and self.__nz != None:
            self.__check_if_same_len((self.__ny, 'ny'), (self.__nz, 'nz'))
        #}}}

        #{{{Check (zperiod), (zmin, zmax) and (dz) is not set simultaneously
        if (self.__zperiod != None and\
           (self.__zmin != None or self.__zmax != None)):
            self.__errors.append("TypeError")
            message = "zperiod and zmin or zmax cannot be set simulatneously."
            raise TypeError(message)
        elif (self.__dz != None and\
             (self.__zmin != None or self.__zmax != None)):
            self.__errors.append("TypeError")
            message = "dz and zmin or zmax cannot be set simulatneously."
            raise TypeError(message)
        elif (self.__zperiod != None and self.__dz):
            self.__errors.append("TypeError")
            message = "dz and zperiod cannot be set simulatneously."
            raise TypeError(message)
        #}}}

        #{{{Check that dz is not set
        # dz is currently set throught zmin and zmax
        if self.__dz != None:
            self.__errors.append("TypeError")
            message  = "dz can currently just be set through zmin and zmax\n"
            message += "dz = 2*pi*(zmax-zmin)/(MZ-1)"
            raise TypeError (message)
        #}}}

        #{{{Check that dx, dy and dz are of the same length
        if self.__dx != None and self.__dy != None:
            self.__check_if_same_len((self.__dx, 'dx'), (self.__dy, 'dy'))
        if self.__dx != None and self.__dz != None:
            self.__check_if_same_len((self.__dx, 'dx'), (self.__dz, 'dz'))
        if self.__dy != None and self.__dz != None:
            self.__check_if_same_len((self.__dy, 'dy'), (self.__dz, 'dz'))
        #}}}

        #{{{Check that (dx, nx), (dy, ny) and (dz,nz) are of the same length
        if self.__dx != None and self.__nx != None:
            self.__check_if_same_len((self.__dx, 'dx'), (self.__nx, 'nx'))
        if self.__dy != None and self.__ny != None:
            self.__check_if_same_len((self.__dy, 'dy'), (self.__ny, 'ny'))
        if self.__nz != None and self.__dz != None:
            self.__check_if_same_len((self.__dz, 'dz'), (self.__nz, 'nz'))
        #}}}

        #{{{ Check that timestep and nout have the same len
        if self.__timestep != None and self.__nout != None:
            self.__check_if_same_len((self.__timestep, 'timestep'),\
                                   (self.__nout, 'nout'))
        #}}}

        #{{{Check that additional and series_add are on the correct form
        self.__error_check_additional((self.__additional, 'additional'))
        self.__error_check_additional((self.__series_add, 'series_add'))
        #}}}

        #{{{Check that self.__series_add[:][2] have the same length
        if self.__series_add != None:
            # Collect all second indices
            second_indices = [elems[2] for elems in self.__series_add]
            # Find the length of the second indices
            lengths = [len(elem) for elem in second_indices\
                       if (type(elem)!=str and type(elem)!=dict)]
            # Check if any string or dicts were given
            if len(second_indices) != len(lengths):
                message  = "self.__series_add is on the wrong form.\n"
                message += "self.__series_add should be on the form\n"
                message += "self.__series_add=\ \n"
                message +=\
                        "     [(section1, name1, [value1-1, value1-2,... ,"+\
                        "value1-n]),\ \n"
                message +=\
                        "      (section2, name2, [value2-1, value2-2,... ,"+\
                        "value2-n]),\ \n"
                message +=\
                        "       ...])\n"
                self.__errors.append("TypeError")
                raise TypeError(message)

            # Check that the length of the second indices are the same
            # L.count(value) -> integer -- return number of occurrences
            # of value
            # stackoverflow.com/questions/3844801/check-if-all-elements-in-a-list-are-identical
            if not(lengths.count(lengths[0]) == len(lengths)):
                message = "The length of the second index of the elements"+\
                          " of self.__series_add must be the same"
                self.__errors.append(message)
                raise TypeError(message)
        #}}}

        #{{{Check mms, symGlobX, symGlobY, cpy_src/grid, allow_size_mod is bool
        check_if_bool = [\
            (self.__mms                    , 'mms')                    ,\
            (self.__symGlobX               , 'symGlobX')               ,\
            (self.__symGlobY               , 'symGlobY')               ,\
            (self.__cpy_source             , 'cpy_source')             ,\
            (self.__cpy_grid               , 'cpy_grid')               ,\
            (self.__allow_size_modification, 'allow_size_modification') \
            ]

        self.__check_for_correct_type(var = check_if_bool,\
                                    the_type = bool)
        #}}}

        #{{{Check grid_file == None if cpy_grid==True
        if (self.__grid_file == None) and (self.__cpy_grid == True):
            # Raise error
            self.__errors.append("TypeError")
            message = "Cannot copy the grid files if none exists in "+\
                      " self.__grid_file"
            raise TypeError(message)
        #}}}

        #{{{Check that zmin and zmax has the same length
        if (self.__zmin != None) and (self.__zmax != None):
            self.__check_if_same_len((self.__zmin, 'zmin'),\
                                     (self.__zmax, 'zmax'))

        #}}}
#}}}
#}}}

#{{{Functions called by __check_for_basic_instance_error
    #{{{
    def __error_check_additional(self, input_member):
        #{{{docstring
        """
        Checks that the input_member is on the following form:

        input_member = [(section1, name1, [value1-1, value1-2, ...]),
                        (section2, name2, [value2-1, value2-2, ...]),
                        ...]

        Input:
        Either self.__additional or self.__series_add
        input_member[0]     -   the input data
        input_member[1]     -   the name of the input data
        """
        #}}}

        # If input_member is set
        if input_member[0] != None:
            # Set a success variable that will fail if anything goes
            # wrong
            success = True
            # Loop through all elements in input_member
            for elem in input_member[0]:
                # Check if self.__addition is iterable, but not a string
                # or dict
                if (hasattr(elem, "__iter__")) and\
                   (type(elem) != str) and\
                   (type(elem) != dict):
                    if type(elem[0]) == str:
                        # Check that the second element (the name) is a
                        # string
                        if type(elem[1]) != str:
                            success = False
                        # If more than three elements are given
                        if len(elem) != 3:
                            success = False
                    # elem[0] is not a string
                    else:
                        success = False
                # elem is not iterable or is a dict or a string
                else:
                    success = False
            if not(success):
                message  = "self.__"+input_member[1]+" is on the wrong form.\n"
                message += "self.__"+input_member[1]+" should be on the form\n"
                message += "self.__"+input_member[1]+"=\ \n"
                message +=\
                        "     [(section1, name1, [value1-1, value1-2,...]),\ \n"
                message +=\
                        "      (section2, name2, [value2-1, value2-2,...]),\ \n"
                message +=\
                        "       ...])\n"
                self.__errors.append("TypeError")
                raise TypeError(message)
        #}}}
#}}}

#{{{ Functions called by the run function
#{{{__error_check_for_run_input
    def __error_check_for_run_input(self,\
                                    remove_old,\
                                    post_processing_function,\
                                    post_process_after_every_run):
        """Check if there are any type errors in input for the run
        function"""

        #{{{Check if remove_old is of the correct type
        check_if_bool = [\
            (remove_old , 'remove_old'),\
            ]

        self.__check_for_correct_type(var = check_if_bool,\
                                    the_type = bool)
        #}}}

        #{{{Check if remove_old and restart is set on the same time
        if remove_old == True and self.__restart != None:
            self.__errors.append("RuntimeError")
            raise RuntimeError("You should not remove old data if you"\
                               " want a restart run")
        #}}}

        #{{{Check that the post_processing_function is a fuction
        if (post_processing_function != None) and\
           (not(hasattr(post_processing_function, '__call__'))):
            self.__errors.append("RuntimeError")
            message = "post_process_after_every_run must be a"+\
                      " function"
            raise RuntimeError(message)
        #}}}

        #{{{Check that the post_process_after_every_run is not set alone
        if (post_process_after_every_run != None) and\
           (type(post_processing_function) == None):
            self.__errors.append("RuntimeError")
            message = "post_process_after_every_run can only be set if"+\
                      " post_processing_function is given"
            raise RuntimeError(message)
        #}}}

        #{{{Check that the post_process_after_every_run is a boolean
        if (post_process_after_every_run != None) and\
           (type(post_process_after_every_run) != bool):
            self.__errors.append("RuntimeError")
            message = "post_process_after_every_run must be set to"+\
                      " a boolean when set"
            raise RuntimeError(message)
        #}}}
#}}}

#{{{__create_run_log
    def __create_run_log(self):
        """Makes a run_log file if it doesn't exists"""

        # Checks if run_log exists
        self.__run_log = self.__directory + "/run_log.txt"
        if os.path.isfile(self.__run_log) == False:
            # The header
            header = ['start_time', 'run_type', 'run_no', 'run_time_H:M:S', 'dump_folder']
            header_format = '{:<19}   {:<9}   {:<6}   {:<17}   {:<}'
            # Create the log file, and print the header
            with open(self.__run_log , "w") as f:
                f.write(header_format.format(*header) + '\n')

        # Preparation of the run
        print("\nRunning with inputs from '" + self.__directory + "'")
#}}}

#{{{__get_correct_domain_split
    def __get_correct_domain_split(self):
        """Checks that the grid can be split in the correct number of
        processors. If not, vary the number of points until value is found."""

        # Flag which is True when a warning should be produce
        produce_warning = False

        # First we check if self.__MXG is given
        if self.__MXG == None:
            # We need to find MXG in BOUT.inp
            local_MXG = self.__find_variable_in_BOUT_inp("MXG")

            # Check that one unique variable was found
            if len(local_MXG) > 1:
                self.__errors.append("RuntimeError")
                message =  "Several matches was found when searching for "
                message += "MXG.\n"
                message += "Please check " + self.__directory + "/BOUT.inp"
                raise RuntimeError(message)

            # Pick the only element
            local_MXG = int(local_MXG[0])
        else:
            local_MXG = self.__MXG

        print("\nChecking the grid split for the meshes\n")
        for size_nr in range(len(self.__nx)):
            print("Checking nx=" + str(self.__nx[size_nr]) +\
                  " and ny=" + str(self.__ny[size_nr]))
            split_found = False
            add_number = 1
            while split_found == False:
                # The same check as below is performed internally in
                # BOUT++ (see boutmesh.cxx)
                for i in range(1, self.__nproc+1, 1):
                    MX = self.__nx[size_nr] - 2*local_MXG
                    if (self.__nproc % i == 0) and \
                       (MX % i == 0) and \
                       (self.__ny[size_nr] % (self.__nproc/i) == 0):
                        # If the test passes
                        split_found = True

                # If the value tried is not a good value
                if split_found == False:
                    # If modification is allowd
                    if self.__allow_size_modification and self.__grid_file == None:
                        # Produce a warning
                        produce_warning = True
                        self.__nx[size_nr] += add_number
                        self.__ny[size_nr] += add_number
                        print("Mismatch, trying "+ str(self.__nx[size_nr]) +\
                              "*" + str(self.__ny[size_nr]))
                        add_number = (-1)**(abs(add_number))\
                                     *(abs(add_number) + 1)
                    else:
                        # If the split fails and the a grid file is given
                        if self.__grid_file != None:
                            self.__errors.append("RuntimeError")
                            message = "The grid can not be split using the"+\
                                      " current number of nproc"
                            raise RuntimeError(message)
                        # If the split fails and no grid file is given
                        else:
                            self.__errors.append("RuntimeError")
                            message  = "The grid can not be split using the"+\
                                       " current number of nproc.\n"
                            message += "Setting allow_size_modification=True"+\
                                       " will allow modification of the grid"+\
                                       " so that it can be split with the"+\
                                       " current number of nproc"
                            raise RuntimeError(message)
            # When the good value is found
            print("Sucessfully found the following good values for the mesh:")
            print("nx=" + str(self.__nx[size_nr]) +\
                  " ny=" + str(self.__ny[size_nr]) +\
                  "\n")

            # Make the warning if produced
            if produce_warning:
                message = "The mesh was changed to allow the split given by nproc"
                self.__warning_printer(message)
                self.__warnings.append(message)
#}}}

#{{{__get_possibilities
    def __get_possibilities(self):
        """ Returns the list of the possibilities. In get_combinations
        the elements of this list is going to be put together to a list
        of strings which will be used when making a run."""

        #{{{Set combination of nx, ny and nz (if not set in gridfile)
        # Appendable list
        spatial_grid_possibilities = []
        if (self.__grid_file == None):
            # Dictionary where
            # - the first element is the variable itself
            # - the second element is the section of the variable
            # - the third element is an appendable list
            spatial_grid_str = {\
                                'nx'     :[self.__nx,      'mesh:', []],\
                                'ny'     :[self.__ny,      'mesh:', []],\
                                'nz'     :[self.__nz,      'mesh:', []],\
                                'dx'     :[self.__dx,      'mesh:', []],\
                                'dy'     :[self.__dy,      'mesh:', []],\
                                'dz'     :[self.__dz,      'mesh:', []],\
                                'zperiod':[self.__zperiod, '',      []],\
                                'zmin'   :[self.__zmin,    '',      []],\
                                'zmax'   :[self.__zmax,    '',      []],\
                               }
            # Store the keys as an own variable
            keys = list(spatial_grid_str.keys())
            # Append the different dimension to the list of strings
            for key in keys:
                # If the variable is not empty
                if spatial_grid_str[key][0] != None:
                    # Fill the appendable list with the elements from
                    # the variable
                    for elem in spatial_grid_str[key][0]:
                        spatial_grid_str[key][2].append(\
                            spatial_grid_str[key][1] + key + '=' + str(elem)\
                                )

            # The goal is to combine the these strings to one string
            # Find the largest length
            lengths = [len(spatial_grid_str[key][2]) for key in keys]
            max_len = np.max(lengths)
            # Make the strings the same length
            for key in keys:
                # We do this by filling it with empty strings
                while len(spatial_grid_str[key][2]) <= max_len:
                    spatial_grid_str[key][2].append('')

            # Append this to the spatial grid possibilities as a string
            for number in range(max_len):
                # Make a list
                current_grid = [spatial_grid_str[key][2][number] for key in\
                                keys]
                # Join the strings in the list and append
                spatial_grid_possibilities.append(' '.join(current_grid))
        #}}}

        #{{{Set the combination of timestep and nout if != None
        # Appendable lists
        temporal_grid_possibilities = []
        timestep_str = []
        nout_str     = []
        # Append the different time options to the list of strings
        if self.__timestep != None:
            for timestep in self.__timestep:
                timestep_str.append('timestep=' + str(timestep))
        if self.__nout != None:
            for nout in self.__nout:
                nout_str.append('nout=' + str(nout))
        # Combine the strings to one string
        # Find the largest length
        max_len = np.max([len(timestep_str), len(nout_str)])
        # Make the strings the same length
        if len(timestep_str) < max_len:
            timestep_str.append('')
        if len(nout_str) < max_len:
            nout_str.append('')
        # Append the temporal grid possibilities as a string
        for number in range(max_len):
            # Make a list
            current_times = [timestep_str[number],\
                             nout_str[number]\
                            ]
            # Join the strings in the list and append
            temporal_grid_possibilities.append(' '.join(current_times))
        #}}}

        #{{{Set the combination of the series_add option if != None
        # Appendable list
        series_add_possibilities = []
        if self.__series_add != None:
            # Dictionary to handle the data, where the key is going to
            # be the element number in self.__series_add, and the values
            # are going to be the sub dictionary defined below
            all_info = {}
            # Loop through all elements and fill the dictionary
            for nr, elem in enumerate(self.__series_add):
                # Put in the sub dictionary
                all_info[nr] = {'values':None,\
                                'section_and_var':None,\
                                'sec_var_vals':[]}
                # Fill the values
                all_info[nr]['values'] = elem[2]
                # Fill the section and variable key
                all_info[nr]['section_and_var'] = elem[0] + ":" + elem[1] + "="
                # Fill in the combinations
                for val in all_info[nr]['values']:
                    all_info[nr]['sec_var_vals'].append(\
                        all_info[nr]['section_and_var'] + str(val)\
                            )

            # Make an appendable list
            all_sec_var_vals = []
            for key in all_info.keys():
                all_sec_var_vals.append(all_info[key]['sec_var_vals'])

            # Zip the sec_var_vals together (* unpacks), join them with
            # a space, and append them to series_add_possibilities
            for one_possibility in zip(*all_sec_var_vals):
                series_add_possibilities.append(' '.join(one_possibility))
        #}}}

        #{{{Put non-iterable variables into a list if they are not set to None
        # This makes the memberdata iterable, and useable in
        # generate_possibilities
        if self.__MXG != None:
            self.__MXG = [self.__MXG]
        if self.__MYG != None:
            self.__MYG = [self.__MYG]
        if self.__mms != None:
            self.__mms = [self.__mms]
        if self.__symGlobX != None:
            self.__symGlobX = [self.__symGlobX]
        if self.__symGlobY != None:
            self.__symGlobY = [self.__symGlobY]
        #}}}

        #{{{List of tuple of varibles to generate possibilities from
        tuple_of_variables = [\
            (self.__solver,     "solver", "type")             ,\
            (self.__mms,        "solver", "mms")              ,\
            (self.__atol,       "solver", "atol")             ,\
            (self.__rtol,       "solver", "rtol")             ,\
            (self.__mxstep,     "solver", "mxstep")           ,\
            (self.__MXG,        "",       "MXG")              ,\
            (self.__MYG,        "",       "MYG")              ,\
            (self.__grid_file,  "",       "grid")             ,\
            (self.__ddx_first,  "ddx",    "first")            ,\
            (self.__ddx_second, "ddx",    "second")           ,\
            (self.__ddx_upwind, "ddx",    "upwind")           ,\
            (self.__ddx_flux,   "ddx",    "flux")             ,\
            (self.__ddy_first,  "ddy",    "first")            ,\
            (self.__ddy_second, "ddy",    "second")           ,\
            (self.__ddy_upwind, "ddy",    "upwind")           ,\
            (self.__ddy_flux,   "ddy",    "flux")             ,\
            (self.__ddz_first,  "ddz",    "first")            ,\
            (self.__ddz_second, "ddz",    "second")           ,\
            (self.__ddz_upwind, "ddz",    "upwind")           ,\
            (self.__ddz_flux,   "ddz",    "flux")             ,\
            (self.__ixseps1,    "mesh",   "ixseps1")          ,\
            (self.__ixseps2,    "mesh",   "ixseps2")          ,\
            (self.__jyseps1_1,  "mesh",   "jyseps1_1")        ,\
            (self.__jyseps1_2,  "mesh",   "jyseps1_2")        ,\
            (self.__jyseps2_1,  "mesh",   "jyseps2_1")        ,\
            (self.__jyseps2_2,  "mesh",   "jyseps2_2")        ,\
            (self.__symGlobX,   "mesh",   "symmetricGlobalX") ,\
            (self.__symGlobY,   "mesh",   "symmetricGlobalY")  \
            ]
        #}}}

        #{{{Append the additional option to tuple of variables if set
        if self.__additional != None:
            for additional in self.__additional:
                # If the last element of additional is not iterable we need
                # put them into a list to make them iterable (in order to
                # use them in generate_possibilities)
                if (not(hasattr(additional[2], "__iter__"))) or\
                   (type(additional[2]) == str):
                    # We have to specify the whole additional, as this can
                    # be given as a tuple, and tuples does not support item
                    # assignments
                    additional = (additional[0],\
                                  additional[1],\
                                  [additional[2]])
                # Append the additional to tuple of variables
                tuple_of_variables.append(\
                                (additional[2],\
                                additional[0],\
                                additional[1])\
                                )
        #}}}

        #{{{List of the possibilities of the variables
        # Start out with the already generated
        # spatial_grid_possibilities and temporal_grid_possibilities
        list_of_possibilities = [spatial_grid_possibilities,\
                                 temporal_grid_possibilities,\
                                 series_add_possibilities]

        # Append the possibilities to the list of possibilities
        for var in tuple_of_variables:
            list_of_possibilities.append(\
                    self.__generate_possibilities(var[0], var[1], var[2])\
                    )
        #}}}

        # Return the list_of possibilities
        return list_of_possibilities
#}}}

#{{{__get_combinations
    def __get_combinations(self, input_list):
        """The input_list is a list with lists as element.
        Returns a list of all combinations between the elements of the
        input_list."""

        # Remove empty elements in input_list in order for
        # itertools.product to work
        input_list = [elem for elem in input_list if elem != []]

        # If we would like to sort the input list (choose which varible
        # to be the fastest varying)
        if self.__sort_by != None:
            # Swap the list corresponding to the sort_by statement so
            # that that list will be the last. The itertools.product
            # will then make that list the fastest varying in the list
            input_list = self.__get_swapped_input_list()
        else:
            # Initialize this member data to None
            self.__len_group = None

        # The last element in the input_list will be the fastest varying
        # element
        all_combinations_as_tuple = list(itertools.product(*input_list))

        # all_combination_as_tuple is a list with tuples as elements
        # We would like to combind the elements in these tuples to one
        # string
        # Make an appendable list
        all_combinations_as_strings = []

        # Loop over the elements in the list containing tuples
        for a_tuple in all_combinations_as_tuple:
            # Join the elements in a tuple and store it
            all_combinations_as_strings.append(' '.join(a_tuple))

        return all_combinations_as_strings
#}}}

#{{{__print_run_or_submit
    def __print_run_or_submit(self):
        """Prints 'Now running'"""
        print("\nNow running:")
#}}}

#{{{__prepare_dmp_folder
    def __prepare_dmp_folder(self, combination):
        """Set the folder to dump data in based on the input from the
        combination. Copy the input file to the final folder. Copy the
        source files to the final folder is cpy_source is True."""
        # Obtain folder names
        folder_name = self.__get_folder_name(combination)
        self.__dmp_folder = os.path.join(self.__directory, folder_name)
        # If the last character is '/', then remove it
        if self.__dmp_folder[-1] == '/':
            self.__dmp_folder = self.__dmp_folder[:-1]

        # Create folder if it doesn't exists
        self.__create_folder(self.__dmp_folder)
        # If self.__dmp_folder contains anything other than
        # self.__directory
        if self.__dmp_folder != self.__directory:
            # Copy the input file into this folder
            command = 'cp ' + self.__directory + '/BOUT.inp ' +\
                      self.__dmp_folder + '/'
            shell(command)

        # Copy the source files if cpy_source is True
        if self.__cpy_source:
            # This will copy all C++ files to the dmp_folder
            cpp_extension= ['.cc', '.cpp', '.cxx', '.C', '.c++',\
                            '.h',  '.hpp', '.hxx', '.h++']
            # Copy for all files in the extension
            for extension in cpp_extension:
                file_names = glob.glob('*' + extension)
                for a_file in file_names:
                    command = 'cp ' + a_file + ' ' + self.__dmp_folder + '/'
                    shell(command)
#}}}

#{{{__remove_data
    def __remove_data(self):
        """Removes *.nc and *.log files from the dump directory"""

        print("Removing old data")
        # Make the command
        command = "rm -f ./" + self.__dmp_folder +\
                  "/*.nc ./" + self.__dmp_folder +\
                  "/*.log.*"
        # Execute the command
        shell(command)
#}}}

#{{{__check_if_run_already_performed
    def __check_if_run_already_performed(self):
        """Checks if the run has been run previously"""

        dmp_files = glob.glob(self.__dmp_folder + '/BOUT.dmp.*')
        # If no BOUT.inp files are found or if self.__restart is not set
        # (meaning that the run will be done even if files are found)
        if len(dmp_files) != 0 and self.__restart == None:
            print('Skipping the run as *.dmp.* files was found in '\
                  + self.__dmp_folder)
            print('To overwrite old files, run with'+\
                  ' self.execute_runs(remove_old=True)\n')
            return False
        else:
            return True
#}}}

#{{{__call_post_processing_function
    def __call_post_processing_function(\
                    self           ,\
                    function = None,\
                    folders  = None,\
                    **kwargs):
        """Function which calls the post_processing_function"""

        function(folders, **kwargs)

#}}}

#{{{__post_run
    def __post_run(self, **kwarg):
        """In basic_runner this is a virtual function"""

        return
#}}}
#}}}

#{{{Function called by __set_program_name
#{{{__run_make
    def __run_make(self):
        """Makes the .cxx program, saves the make.log and make.err"""

        print("Making the .cxx program\n")
        command = "make"
        status, output = shell(command, pipe=True)
        print(output)
        # Check if any errors occured
        if status != 0:
            self.__errors.append("RuntimeError")
            raise RuntimeError("Error encountered during make.")
#}}}
#}}}

#{{{ Functions called by the basic_error_checker
#{{{__check_for_correct_type
    def __check_for_correct_type(self,\
                                 var            = None,\
                                 the_type       = None,\
                                 allow_iterable = None):
        """Checks if a varible has the correct type

        Input:
        var            - a tuple consisting of
                         var[0] - the variable (a data member)
                         var[1] - the name of the varible given as a string
        the_type       - the data type
        allow_iterable - if an iterable with the element as type is
                         allowed
        """

        # Set a variable which is False if the test fails
        success = True
        for cur_var in var:
            # There is an option that the variable could be set to None,
            # and that the default value from BOUT.inp will be used
            if cur_var[0] != None:
                # Check for the correct type
                if isinstance(cur_var[0], the_type) == False:
                    # Check if it is an iterable if iterables are
                    # allowed
                    if allow_iterable and\
                       hasattr(cur_var[0], "__iter__") and\
                       type(cur_var[0]) != dict:
                        for elem in cur_var[0]:
                            # Check for the correct type
                            if isinstance(elem, the_type) == False:
                                success = False
                    else:
                        # Neither correct type, nor iterable
                        success = False
                if not(success):
                    message  = cur_var[1] + " is of wrong type\n"+\
                               cur_var[1] + " must be " + the_type.__name__
                    if allow_iterable:
                        # If iterable is allowed, then add this
                        message += " or an iterable with " + the_type.__name__
                        message += " as elements."
                    self.__errors.append("TypeError")
                    raise TypeError(message)
#}}}

#{{{__check_if_set_correctly
    def __check_if_set_correctly(self,\
                                 var           = None,\
                                 possibilities = None):
        """Check if a variable is set to a possible variable. Called by
        the error checkers"""

        # Set a variable which is False if the test fails
        success = True

        # Due to the check done in check_for_correct_type: If the
        # variable is not a string it will be an iterable
        if type(var[0]) != str:
            for elem in var[0]:
                # Check if the element is contained in the possibilities
                if not(elem in possibilities):
                    success = False
        else:
            # The variable was a string
            if not(var[0] in possibilities):
                success = False

        if not(success):
            message = var[1] + " was not set to a possible option.\n"+\
                      "The possibilities are \n" + "\n".join(possibilities)
            self.__errors.append("TypeError")
            raise TypeError(message)
#}}}

#{{{__check_if_same_len
    def __check_if_same_len(self, object1 = None, object2 = None):
        """Checks if object1 and object2 has the same length

        Input:
        object1 - a tuple of the object [0] and its name [1]
        object2 - a tuple an object [0] different than object1 together with
                  its name [1]
        """

        try:
            len_dim1 = len(object1[0])
        # If object1 does not have length
        except TypeError:
            len_dim1 = 1
        try:
            len_dim2 = len(object2[0])
        # If object2 does not have length
        except TypeError:
            len_dim2 = 1

        if len_dim1 != len_dim2:
            message = object1[1] + " and " + object2[1] + " must have the same"
            message += " length when specified"
            self.__errors.append("RuntimeError")
            raise RuntimeError (message)
#}}}
#}}}

#{{{Function called by __prepare_dmp_folder
#{{{__get_folder_name
    def __get_folder_name(self, combination):
        """Returning the folder name where the data will be stored.

        If all options are given the folder structure should be on the
        form solver/method/nout_timestep/mesh/additional/grid"""

        # Combination is one of the combination of the data members
        # which is used as the command line arguments in the run
        combination = combination.split()

        #{{{Append from eventual grid file
        # If there is a grid file, we will extract the values from the
        # file, and put it into this local combination variable, so that
        # a proper dmp folder can be made on basis on the variables
        # A flag to see whether or not the grid file was found
        grid_file_found = False
        # Check if grid is in element, and extract its path
        for elem in combination:
            if 'grid' in elem:
                cur_grid = elem.replace('grid=','')
                grid_file_found = True

        # If the grid file is found, open it
        if grid_file_found:
            # Open (and automatically close) the grid files
            f = DataFile(cur_grid)
            # Search for mesh types in the grid file
            mesh_types = [\
                          ("mesh:", "nx"       ) ,\
                          ("mesh:", "ny"       ) ,\
                          ("mesh:", "nz"       ) ,\
                          ("mesh:", "zperiod"  ) ,\
                          ("mesh:", "zmin"     ) ,\
                          ("mesh:", "zmax"     ) ,\
                          ("mesh:", "dx"       ) ,\
                          ("mesh:", "dy"       ) ,\
                          ("mesh:", "dz"       ) ,\
                          ("mesh:", "ixseps1"  ) ,\
                          ("mesh:", "ixseps2"  ) ,\
                          ("mesh:", "jyseps1_1") ,\
                          ("mesh:", "jyseps1_2") ,\
                          ("mesh:", "jyseps2_1") ,\
                          ("mesh:", "jyseps2_2") ,\
                          ("",      "MXG")       ,\
                          ("",      "MYG")       ,\
                          ]
            for mesh_type in mesh_types:
                grid_variable = f.read(mesh_type[1])
                # If the variable is found
                if grid_variable != None:
                    # Append it to the combinations list
                    combination.append(mesh_type[0] +\
                                       mesh_type[1] +\
                                       "=" + str(grid_variable))
        #}}}

        # Make lists for the folder-type, so that we can append the
        # elements in the combination folders if it is found
        solver        = []
        method        = []
        nout_timestep = []
        mesh          = []
        additional    = []
        grid_file     = []

        # We will loop over the names describing the methods used
        # Possible directional derivatives
        dir_derivatives = ['ddx', 'ddy', 'ddz']

        # Check trough all the elements of combination
        for elem in combination:

            # If 'solver' is in the element
            if 'solver' in elem:
                # Remove 'solver:' and append it to the mesh folder
                cur_solver = elem.replace('solver:','')
                cur_solver = cur_solver.replace('=','_')
                # Append it to the solver folder
                solver.append(cur_solver)

            # If nout or timestep is in the element
            elif ('nout' in elem) or\
                 ('timestep' in elem):
                # Remove '=', and append it to the
                # nout_timestep folder
                nout_timestep.append(elem.replace('=','_'))

            # If mesh, MXG or MYG is in the combination
            elif ('mesh' in elem) or\
                 ('MXG' in elem) or\
                 ('MYG' in elem) or\
                 ('zperiod' in elem) or\
                 ('zmin' in elem) or\
                 ('zmax' in elem) or\
                 (('dx' in elem) and not('ddx' in elem)) or\
                 (('dy' in elem) and not('ddy' in elem)) or\
                 (('dz' in elem) and not('ddz' in elem)):
                # Remove 'mesh:', and append it to the mesh folder
                cur_mesh = elem.replace('mesh:','')
                cur_mesh = cur_mesh.replace('=','_')
                mesh.append(cur_mesh)

            # If a grid file is in the combination
            elif ('grid' in elem):
                # Remove .nc and =
                cur_grid = elem.replace('.nc','')
                cur_grid = cur_grid.replace('=','_')
                grid_file.append(cur_grid)

            # If the element is none of the above
            else:
                # It could either be a dir derivative
                # Set a flag to state if any of the dir derivative was
                # found in the combination
                dir_derivative_set = False
                # If any of the methods are in combiniation
                for dir_derivative in dir_derivatives:
                    if dir_derivative in elem:
                        # Remove ':', and append it to the
                        # method folder
                        cur_method = elem.replace(':','_')
                        cur_method = cur_method.replace('=','_')
                        method.append(cur_method)
                        dir_derivative_set = True

                # If the dir_derivative_set was not set, the only
                # possibility left is that the element is an
                # 'additional' option
                if not(dir_derivative_set):
                    # Replace ':' and '=' and append it to the
                    # additional folder
                    cur_additional = elem.replace(':','_')
                    cur_additional = cur_additional.replace('=','_')
                    additional.append(cur_additional)

        # We sort the elements in the various folders alphabethically,
        # to ensure that the naming convention is always the same, no
        # matter how the full combination string looks like
        # Sort alphabethically
        solver.sort()
        #{{{ Manual sort solver
        # We want 'type' to be first, and 'atol' and 'rtol' to be last
        sort_these = [\
                      ('type',0) ,\
                      ('atol',-1),\
                      ('rtol',-1) \
                     ]
        # Loop through everything we want to sort
        for sort_this in sort_these:
            # Flag to check if found
            found_string = False
            for elem_nr, elem in enumerate(solver):
                if sort_this[0] in elem:
                    swap_nr = elem_nr
                    # Set the flag that the string is found
                    found_string = True
            # If type was found
            if found_string != False:
                # Swap the elements in the solver
                solver[sort_this[1]], solver[swap_nr] =\
                        solver[swap_nr], solver[sort_this[1]]
        #}}}
        method.sort()
        nout_timestep.sort()
        mesh.sort()
        additional.sort()
        grid_file.sort()

        # Combine the elements in the various folders
        solver              = ['_'.join(solver)]
        method              = ['_'.join(method)]
        nout_timestep       = ['_'.join(nout_timestep)]
        mesh                = ['_'.join(mesh)]
        additional          = ['_'.join(additional)]
        grid_file           = ['_'.join(grid_file)]

        # Put all the folders into the combination_folder
        combination_folder = [\
                              solver       ,\
                              method       ,\
                              nout_timestep,\
                              mesh         ,\
                              additional   ,\
                              grid_file     \
                             ]
        # We access the zeroth element (if given) as the folders are
        # given as a list
        combination_folder = [folder[0] for folder in combination_folder\
                              if (folder != []) and (folder !=[''])]

        # Make the combination folder as a string
        combination_folder = '/'.join(combination_folder)

        return combination_folder
#}}}

#{{{__create_folder
    def __create_folder(self, folder):
        """Creates a folder if it doesn't exists"""

        if not os.path.exists(folder):
            os.makedirs(folder)
            print(folder + " created\n")
#}}}
#}}}

#{{{Function called by __run_driver
#{{{__single_run
    def __single_run(self, combination):
        """Makes a single MPIRUN of the program"""

        # Get the command to be used
        command = self.__get_command_to_run(combination)

        # Time how long the time took
        tic = timeit.default_timer()

        # Launch the command
        status, out = launch(command,\
                             runcmd = self.__MPIRUN,\
                             nproc = self.__nproc,\
                             pipe = True,\
                             verbose = True)

        # If the run returns an exit code other than 0
        if status != 0:
            message = "! An error occured. Printing the output to stdout !"
            print("\n" + "!"*len(message))
            print(message)
            print("!"*len(message) + "\n")
            print(out)
            self.__errors.append("RuntimeError")
            message =  "An error occured the run."
            message += " Please see the output above for details."
            raise RuntimeError(message)

        # Estimate elapsed time
        toc = timeit.default_timer()
        elapsed_time = toc - tic

        return out, elapsed_time
#}}}

#{{{__append_run_log
    def __append_run_log(self, start, run_no, run_time):
        """Appends the run_log"""

        # Convert seconds to H:M:S
        run_time = str(datetime.timedelta(seconds=run_time))

        start_time = (str(start.year) + '-' + str(start.month) + '-' +\
                      str(start.day) + '.' + str(start.hour) + ":" +\
                      str(start.minute) + ':' + str(start.second))

        # If the run is restarted with initial values from the last run
        if self.__restart:
            dmp_line = self.__dmp_folder + '-restart-'+self.__restart
        else:
            dmp_line = self.__dmp_folder

        # Line to write
        line = [start_time, self.__run_type, run_no, run_time, dmp_line]
        # Opens for appending
        log_format = '{:<19}   {:^9}   {:^6}   {:<17}   {:<}'
        with open(self.__run_log , "a") as f:
            f.write(log_format.format(*line) + '\n')
#}}}
#}}}

#{{{Function called by __get_correct_domain_split
#{{{__find_variable_in_BOUT_inp
    def __find_variable_in_BOUT_inp(self, variable):
        """Find a variable in BOUT.inp using regex.
        The value of the variable is returned"""

        # First, read the file
        with open (self.__directory + "/BOUT.inp", "r") as BOUT_inp_file:
                BOUT_inp = BOUT_inp_file.readlines()

        # Make the string to search for
        # http://www.tutorialspoint.com/python/python_reg_expressions.htm
        # r denotes a raw string
        # First ^ is the start of the string (excludes comments etc.)
        # \s is any whitespace, * repeats any number of times, () encloses the statement
        # = is the literal =
        # \d is digit
        # {n,m} Matches at least n and at most m occurrences of preceding
        # expression.
        # \. is literal . (as only . matches any single character except newline)
        string_to_match = r'^' + variable +\
                          r'(\s*)=(\s*)(\d*)(\.){0,1}(\d*)'

        # Appendable list
        matches = []
        for line in BOUT_inp:
            matchObj = re.match(string_to_match, line, flags=0)
            if matchObj != None:
                matches.append(matchObj.group())

        # Check that a match was found
        if len(matches) == 0:
            self.__errors.append("RuntimeError")
            message =  "No matches was found when searching for "
            message += str(variable) + ".\n"
            message += "Please check " + self.__directory + "/BOUT.inp"
            raise RuntimeError(message)

        # Make an appendable list
        matches_as_floats = []

        for match_as_string in matches:
            # Remove the preceeding numbers
            # re.sub(search, replace, string, max=0)
            # $ matches the beginning of the line
            # . matches any single character except newline
            # * repeats any number of times
            number_as_string = re.sub(r'^(.*)(\s*)=(\s*)', "", match_as_string)

            # Cast the string to a float and append it
            matches_as_floats.append(float(number_as_string))

        return matches_as_floats
#}}}
#}}}

#{{{Function called by __get_possibilities
#{{{__generate_possibilities
    def __generate_possibilities(self, variables=None, section=None, name=None):
        """Generate the list of strings of possibilities"""

        if variables != None:
            # Set the section name correctly
            if section != "":
                section = section + ":"
            else:
                section = ""
            # Set the combination of the varibale
            var_possibilities = []
            # Find the number of different dimensions
            for var in variables:
                var_possibilities.append(section + name + '=' + str(var))
        else:
            var_possibilities = []

        return var_possibilities
#}}}
#}}}

#{{{Functions called by __get_combinations
#{{{__get_swapped_input_list
    def __get_swapped_input_list(self, input_list):
        """Finds the element in the input list, which corresponds to the
        self.__sort_by criterion. The element is swapped with the last
        index, so that itertools.product will make this the fastest
        varying varibale"""

        # We make a sort list containing the string to find in the
        # input_list
        sort_list = []
        # We also make a list over text to write to the error riser if
        # the text in input_list was not found
        text_if_not_found_list = []

        # We loop over the elements in self.__sort_by to find what
        # string we need to be looking for in the elements of the lists
        # in input_list
        for sort_by in self.__sort_by:
            # Find what list in the input_list which contains what we
            # would sort by

            # If we would like to sort by the spatial domain
            if sort_by == 'spatial_domain':
                # nx, ny and nz are all under the section 'mesh'
                find_in_list = 'mesh'
                # Text to be printed if find_in_list is not found
                text_if_not_found = 'neither nx, ny or nz was found'

            # If we would like to sort by the temporal domain
            elif sort_by == 'temporal_domain':
                # If we are sorting by the temporal domain, we can either
                # search for timestep or nout
                if self.__timestep != None:
                    find_in_list = 'timestep'
                elif self.__nout != None:
                    find_in_list = 'nout'
                else:
                    message = "Could not sort by 'temporal_domain' as"+\
                              " neither 'timestep' nor 'nout' is given"
                    raise RuntimeError(message)

            # If we would like to sort by the method
            elif (sort_by == 'ddx_first') or\
                 (sort_by == 'ddx_second') or\
                 (sort_by == 'ddx_upwind') or\
                 (sort_by == 'ddx_flux') or\
                 (sort_by == 'ddy_first') or\
                 (sort_by == 'ddy_second') or\
                 (sort_by == 'ddy_upwind') or\
                 (sort_by == 'ddy_flux') or\
                 (sort_by == 'ddz_first') or\
                 (sort_by == 'ddz_second') or\
                 (sort_by == 'ddz_upwind') or\
                 (sort_by == 'ddz_flux'):
                find_in_list = sort_by.replace('_',':')
                text_if_not_found = sort_by + " was not found"

            # If we would like to sort by the solver
            elif sort_by == 'solver':
                find_in_list = sort_by
                text_if_not_found = sort_by + "was not found"

            # If we would like to sort by additional
            else:
                find_in_list = sort_by
                text_if_not_found = sort_by + " was not given as "+ \
                                    " an additional option."

            # Append what to be found in the input_list, and what text
            # to be printed if it is not found to the respective lists
            sort_list.append(find_in_list)
            text_if_not_found_list.append(text_if_not_found)


        # For all the find_in_list, we would like check if the match
        # can be found in any of the elements in input_list
        for sort_nr, sort_by_txt in enumerate(sort_list):
            # Make a flag to break the outermost loop if find_in_list is
            # found
            break_outer = False
            # Make a flag to state whether or not the sort criteria was
            # found in the elements of input_list
            found = False
            # Loop over the lists in the input_list to find the match
            for elem_nr, elem in enumerate(input_list):
                # Each of the elements in this list is a string
                for string in elem:
                    # Check if fins_in_list is in the string
                    if sort_by_txt in string:
                        # If there is a match, store the element number
                        swap_from_index = elem_nr
                        # Check the length of the element (as this is
                        # the number of times the run is repeated, only
                        # changing the values of sort_by [defining a
                        # group])
                        len_cur_input_list_elem = len(elem)
                        # Break the loop to save time
                        found = True
                        break_outer = True
                        break
                # Break the outer loop if find_in_list_is_found
                if break_outer:
                    break

            # If there was no match, throw an error
            if not(found):
                message  = "Could not sort by " + sort_by_txt + " as "
                message += text_if_not_found_list[sort_nr]
                raise RuntimeError(message)

            # As it is the last index which changes the fastest, we swap the
            # element where the find_in_list was found with the last element
            input_list[swap_from_index], input_list[-(sort_nr + 1)] =\
                    input_list[-(sort_nr + 1)], input_list[swap_from_index]

        # The number of runs in one 'group'
        # Initialize self.__len_group with one as we are going to
        # multiply it with all the elements in len_cur_input_list_elem
        self.__len_group = 1
        for elem in len_cur_input_list_elem:
            self.__len_group *= len_cur_input_list_elem

        return input_list
#}}}
#}}}

#{{{Function called by __single_run
#{{{__get_command_to_run
    def __get_command_to_run(self, combination):
        """ Returns a string of the command which will run the BOUT++
        program"""

        # Creating the arguments
        arg = " -d " + self.__dmp_folder + " " + combination

        # If the run is restarted with initial values from the last run
        if self.__restart != None:
            if self.__restart == 'overwrite':
                arg += ' restart'
            elif self.__restart == 'append':
                arg += ' restart append'
            else:
                self.__errors.append("TypeError")
                raise TypeError ("self.__restart must be set to either"+\
                                 " 'overwrite' or 'append'")

        # Replace excessive spaces with a single space
        arg = ' '.join(arg.split())
        command = "./" + self.__program_name + " " + arg

        return command
#}}}
#}}}

#{{{Functions called from several places in the code
    #{{{__warning_printer
    def __warning_printer(self, message):
        """Function for printing warnings"""

        print('\n'*3 + '*'*37 + 'WARNING' + '*'*36)
        # Makes sure that no more than 80 characters are printed out at
        # the same time
        message_chunks=[]
        for chunk in self.__message_chunker(message):
            rigth_padding = ' '*(76 - len(chunk))
            print('* ' + chunk + rigth_padding + ' *')
        print('*'*80 + '\n'*3)
    #}}}

    #{{{__message_chunker
    def __message_chunker(self, message, chunk=76):
        """Generator used to chop a message so it doesn't exceed some
        width"""

        for start in range(0, len(message), chunk):
            yield message[start:start + chunk]
    #}}}
#}}}
#}}}



#{{{PBS_runner
class PBS_runner(basic_runner):
#{{{docstring
    """Class for running BOUT++.
    Works like the basic_runner, but submits the jobs to a torque queue
    with qsub.

    The link below gives a nice introduction to the PBS system
    http://wiki.ibest.uidaho.edu/index.php/Tutorial:_Submitting_a_job_using_qsub"""
#}}}

# The constructor
#{{{__init__
    def __init__(self,\
                 nodes      = '1',\
                 ppn        = '4',\
                 walltime   = '50:00:00',\
                 mail       = False,\
                 queue      = False,\
                 solver     = False,\
                 nproc      = 1,\
                 methods    = False,\
                 n_points   = False,\
                 directory  = 'data',\
                 nout       = False,\
                 timestep   = False,\
                 MYG        = False,\
                 MXG        = False,\
                 additional = False,\
                 restart    = False,\
                 **kwargs):
        """The values in the constructor determines the torque job
        is submitted."""

        # Note that the constructor accepts additional keyword
        # arguments. This is because the constructor can be called with
        # 'super' from PBS_run_with_plots, which inherits from both
        # basic_qsub_runner and run_with_plots (which takes different
        # arguments as input)

        # Call the constructor of the superclass
        super(basic_qsub_runner, self.__init__(solver     = solver,\
                                                nproc      = nproc,\
                                                methods    = methods,\
                                                n_points   = n_points,\
                                                directory  = directory,\
                                                nout       = nout,\
                                                timestep   = timestep,\
                                                MYG        = MYG,\
                                                MXG        = MXG,\
                                                additional = additional,\
                                                restart    = restart,\
                                                **kwargs))

        self.__nodes      = nodes
        # Processors per node
        self.__ppn        = ppn
        self.__walltime   = walltime
        self.__mail       = mail
        self.__queue      = queue
        self.__run_type   = 'basic_PBS'
        # A string which will be used to write a self.__deleting python
        # script
        self.__python_tmp = ''
        # The jobid returned from the PBS
        self.__PBS_id = None
#}}}

# The run_driver
#{{{run_driver
    def run_driver(self, do_run, combination, run_no):
        """The machinery which actually performs the run"""
        if do_run:
            job_name = self.__single_submit(combination, run_no)
            # Append the job_name to job_status
            self.__run_groups[self.__group_no]['job_status'].append(job_name)
        else:
            self.__run_groups[self.__group_no]['job_status'].append('done')
#}}}

# Functions called directly by the main function
#{{{
#{{{print_run_or_submit
    def print_run_or_submit(self):
        """Prints 'Submitting'"""
        print("\nSubmitting:")
#}}}

#{{{additional_error_check
    def additional_error_check(self, **kwargs):
        """Calls all the error checkers"""
        self.__PBS_error_check(**kwargs)
#}}}

#{{{single_submit
    def single_submit(self, combination, run_no):
        """Single PBS submission"""
        # Get the script (as a string) which is going to be
        # submitted
        job_name, job_string =\
            self.__get_job_string(run_no, combination)

        # The submission
        self.__PBS_id = self.__submit_to_PBS(job_string)
        return job_name
#}}}

#{{{post_run
    def post_run(self):
        """Creates a self.__deleting python scripts which calls
        clean_up_runs.

        If we would not submit this a job, it would have caused a bottle
        neck if the driver running the basic_runner class would iterate
        over several folders."""
        # Creates a folder to put the .log and .err files created by the
        # PBS in
        create_folder(self.__directory + '/PBS_output')

        # Get the start_time
        start_time = self.__get_start_time()

        # The name of the file
        python_name = 'clean_up_'+start_time+'.py'

        # Creating the job string
        job_name = 'clean_up_' + self.__run_type + '_'+ start_time

        # Get the core of the job_string (note that we only need to use
        # one node and one processor for this)
        job_string = self.__create_PBS_core_string(\
            job_name, '1', '1', self.__walltime,\
            folder = self.__directory + '/PBS_output/')
        # We will write a python script which calls the
        # relevant bout_plotter

        # First line of the script
        self.__python_tmp =\
            'import os\n' +\
            'from bout_runners.common_bout_functions import '+\
            'clean_up_runs\n'
        # Call clean_up_runs
        self.__python_tmp +=\
            "clean_up_runs("+\
            str(self.__run_groups) + ","+\
            "'" + str(self.__directory)  + "')\n"
        # When the script has run, it will delete itself
        self.__python_tmp += "os.remove('" + python_name + "')\n"

        # Write the python script
        f = open(python_name, "w")
        f.write(self.__python_tmp)
        f.close()

        # Call the python script in the submission
        job_string += 'python ' + python_name + '\n'
        job_string += 'exit'

        # Submit the job
        print('\nSubmitting a script which waits for the runs to finish')
        self.__submit_to_PBS(job_string, dependent_job = self.__PBS_id)
#}}}
#}}}

# Auxiliary functions
#{{{
#{{{PBS_error_check
    def PBS_error_check(self, **kwargs):
        """Checks for specific PBS errors"""
        variables = [self.__nodes, self.__ppn, self.__walltime, self.__mail]
        for variable in variables:
            if variable == False:
                # Check that the non-optional variables are set
                if variable == self.__nodes or\
                   variable == self.__ppn or\
                   variable == self.__walltime:
                    if variable == self.__nodes:
                        name = 'self.__nodes'
                    elif variable == self.__ppn:
                        name = 'self.__ppn'
                    elif variable == self.__walltime:
                        name = 'self.__walltime'
                    self.__errors.append("TypeError")
                    raise TypeError (name + " cannot be 'False'.")
            if variable != False:
                # Check that the variables are all given as strings
                if type(variable) != str:
                    self.__errors.append("TypeError")
                    raise TypeError ("All non-false data members in"\
                                     " PBS runners must be strings")
                if variable == self.__nodes:
                    try:
                        int(variable)
                    except ValueError:
                        self.__errors.append("ValueError")
                        raise ValueError ("self.__nodes must be given"\
                                          " as a string of an integer")
                elif variable == self.__ppn:
                    try:
                        int(variable)
                    except ValueError:
                        self.__errors.append("ValueError")
                        raise ValueError ("self.__ppn must be given"\
                                          " as a string of an integer")
                elif variable == self.__walltime:
                    message = "self.__walltime must be on the form 'HH:MM:SS'"
                    # Check if it is given on the form HH:MM:SS
                    walltime_list = self.__walltime.split(':')
                    if len(walltime_list) != 3:
                        self.__errors.append("ValueError")
                        raise ValueError (message)
                    for walltime_no, walltime_element in enumerate(walltime_list):
                        try:
                            int(walltime_element)
                        except ValueError:
                            self.__errors.append("ValueError")
                            raise ValueError (message)
                        # Minutes and seconds can max be 60
                        if walltime_no >= 1 and int(walltime_element) > 60:
                            self.__errors.append("ValueError")
                            raise ValueError (message)
                elif variable == self.__mail:
                    if ('@' in variable) == False and\
                    ('.' in variable) == False:
                        self.__errors.append("ValueError")
                        raise ValueError ("self.__mail must be an email"\
                                          "address")
#}}}

#{{{get_job_string
    def get_job_string(self, run_no, combination):
        """Make a string which will act as a shell script when sent to
        PBS."""

        # Find the combination name
        # Split the name to a list
        combination_name = combination.split(' ')
        # Remove whitespaces
        combination_name = [element for element in combination_name\
                            if element != '']
        # Collect the elements
        combination_name = '_'.join(combination_name)
        # Replace bad characters
        combination_name = combination_name.replace(':','')
        combination_name = combination_name.replace('=','-')

        # Name of job
        job_name = combination_name + '_' + self.__directory + '_' + str(run_no)

        command = self.__get_command_to_run( combination )
        command = 'mpirun -np ' + str(self.__nproc) + ' ' + command

        # Print the command
        print(command + '\n')

        # Get the time for start of the submission
        start = datetime.datetime.now()
        start_time = (str(start.year) + '-' + str(start.month) + '-' +\
                      str(start.day) + '.' + str(start.hour) + ":" +\
                      str(start.minute) + ':' + str(start.second))

        # Creating the job string
        job_string = self.__create_PBS_core_string(\
            job_name, self.__nodes, self.__ppn, self.__walltime)

        # Start the timer
        job_string += 'start=`date +%s`\n'
        # Run the bout program
        job_string += command + '\n'
        # end the timer
        job_string += 'end=`date +%s`\n'
        # Find the elapsed time
        job_string += 'time=$((end-start))\n'
        # The string is now in seconds
        # The following procedure will convert it to H:M:S
        job_string += 'h=$((time/3600))\n'
        job_string += 'm=$((($time%3600)/60))\n'
        job_string += 's=$((time%60))\n'
        # Ideally we would check if any process were writing to
        # run_log.txt
        # This could be done with lsof command as described in
        # http://askubuntu.com/questions/14252/how-in-a-script-can-i-determine-if-a-file-is-currently-being-written-to-by-ano
        # However, lsof is not available on all clusters
        job_string += "echo '" +\
                      start_time + " "*4 +\
                      self.__run_type + " "*4 +\
                      str(run_no) + " "*4 +\
                      self.__dmp_folder + " "*4 +\
                      "'$h':'$m':'$s"+\
                      " >> $PBS_O_WORKDIR/" + self.__directory +\
                      "/run_log.txt \n"
        # Exit the qsub
        job_string += 'exit'

        return job_name, job_string
#}}}

#{{{get_start_time
    def get_start_time(self):
        """Returns a string of the current time down to micro precision"""
        # The time is going to be appended to the  job name and python name
        time_now = datetime.datetime.now()
        start_time = str(getattr(time_now, 'hour')) + '-' +\
                     str(getattr(time_now,'minute'))+ '-' +\
                     str(getattr(time_now,'second'))
        # In case the process is really fast, so that more than one job
        # is submitted per second, we add a microsecond in the
        # names for safety
        start_time += '-' + str(getattr(time_now,'microsecond'))
        return start_time
#}}}

#{{{create_PBS_core_string
    def create_PBS_core_string(\
        self, job_name, nodes, ppn, walltime, folder=''):
        """Creates the core of a PBS script as a string"""

        # Shebang line
        job_string = '#!/bin/bash\n'
        # The job name
        job_string += '#PBS -N ' + job_name + '\n'
        job_string += '#PBS -l nodes=' + nodes + ':ppn=' + ppn  + '\n'
        # Wall time, must be in format HOURS:MINUTES:SECONDS
        job_string += '#PBS -l walltime=' + walltime + '\n'
        if self.__queue != False:
            job_string += '#PBS -q ' + self.__queue + '\n'
        job_string += '#PBS -o ' + folder + job_name + '.log' + '\n'
        job_string += '#PBS -e ' + folder + job_name + '.err' + '\n'
        if self.__mail != False:
            job_string += '#PBS -M ' + self.__mail + '\n'
        # #PBS -m abe
        # a=aborted b=begin e=ended
        job_string += '#PBS -m e ' + '\n'
        # cd to the folder you are sending the qsub from
        job_string += 'cd $PBS_O_WORKDIR ' + '\n'

        return job_string
#}}}

#{{{submit_to_PBS
    def submit_to_PBS(self, job_string, dependent_job=None):
        """Saves the job_string as a shell script, submits it and
        deletes it. Returns the output from PBS as a string"""
        # We will use the subprocess.check_output in order to get the
        # jobid number
        # http://stackoverflow.com/questions/2502833/store-output-of-subprocess-popen-call-in-a-string

        # Create the name of the temporary shell script
        # Get the start_time
        start_time = self.__get_start_time()
        script_name = 'tmp_'+start_time+'.sh'

        # Save the string as a script
        with open(script_name, "w") as shell_script:
                shell_script.write(job_string)

        if dependent_job==None:
            output = check_output(["qsub", "./"+script_name])
        else:
            # http://stackoverflow.com/questions/19517923/how-to-wait-for-a-torque-job-array-to-complete
            output = check_output(["qsub", "depend=afterok:"+dependent_job,\
                                    "./" + script_name])
        # Trims the end of the output string
        output = output.strip(' \t\n\r')

        # Delete the shell script
        command = "rm -f "+script_name
        shell(command)

        return output
#}}}
#}}}
#}}}



#{{{if __name__ == '__main__':
if __name__ == '__main__':
    """If bout_runners is run as a script, it will just call the demo
    function"""

    print("\n\nTo find out about the bout_runners, please read the user's"+\
          " manual, or have a look at 'BOUT/examples/bout_runners_example'")
#}}}
