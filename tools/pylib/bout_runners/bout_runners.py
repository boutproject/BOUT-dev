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
__version__ = '0.9beta'
__date__    = '26.06.2015'


# TODO: SHOULD KEEP PLOTTING ROUTINES? BOUTUTILS?

# TODO: Check if you can delete these
from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
from builtins import object

import textwrap
import os
import re
import itertools
import glob
import timeit
import datetime
import math
import six
from numpy import logspace
from numbers import Number
import numpy as np
from subprocess import check_output
from boututils import shell, launch, getmpirun
from boututils.datafile import DataFile
try:
    # TODO: ACTUALLY DELETE THIS?
    from bout_runners.bout_plotters import convergence_plotter,\
                                           solution_plotter,\
                                           solution_and_error_plotter
    # TODO: ACTUALLY DELETE THIS?
    from bout_runners.common_bout_functions import create_folder,\
                                                   find_variable_in_BOUT_inp,\
                                                   warning_printer,\
                                                   check_for_plotters_errors,\
                                                   clean_up_runs,\
                                                   message_chunker
except ImportError:
    print("Could not load bout_runners from the current directory")

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
    Calling self.run() will run your BOUT++ program with all possible
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

# The constructor
#{{{__init__
    def __init__(self,\
                 nproc      = 1,\
                 directory  = 'data',\
                 solver     = None,\
                 nx         = None,\
                 ny         = None,\
                 nz         = None,\
                 grid_file  = None,\
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
                 MXG        = None,\
                 MYG        = None,\
                 additional = None,\
                 restart    = None,\
                 cpy_source = None,\
                 sort_by    = None,\
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
        nx          -    Number of nx in the run (int or iterable)
        ny          -    Number of ny in the run (int or iterable)
        nz          -    Number of nz in the run (int or iterable)
        grid_file   -    The grid file (str or iterable)
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
        MXG         -    The number of guard cells in the x direction
                         (int)
        MYG         -    The number of guard cells in the y direction
                         (int)
        additional  -    Additional option for the run given on the form
                         ('section_name','variable name', values) or as
                         iterable on the same form, where values can be
                         any value or string or an iterable of those
        restart     -    Wheter or not to use the restart files
                         ('overwrite' or 'append')
        cpy_source  -    Wheter or not to copy the source files to the
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
                        or any 'variable_name' from additional

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
        self.__nx         = self.__set_member_data(nx)
        self.__ny         = self.__set_member_data(ny)
        self.__nz         = self.__set_member_data(nz)
        self.__grid_file  = self.__set_member_data(grid_file)
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
        self.__MXG        = MXG
        self.__MYG        = MYG
        self.__additional = additional
        self.__restart    = restart
        self.__cpy_source = cpy_source
        self.__sort_by    = sort_by
        self.__allow_size_modification = allow_size_modification

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

        # Initializing self.__warnings and self.__error
        # self.__warnings will be filled with warnings
        # self.__errors will be filled with errors
        # The warnings and errors will be printed when the destructor is called
        self.__warnings   = []
        self.__errors     = []

        # Set self.__program_name from the *.o file. Make the program if
        # the *.o file is not found
        self.__set_program_name()
        # Find all files with the extension .o

        # Obtain the MPIRUN
        self.__MPIRUN = getmpirun()

        # Data members used internally in this class and its subclasses
        # The dmp_folder is the the folder where the runs are stored
        # It will be set by self.__prepare_dmp_folder
        self.__dmp_folder = None

        #{{{TODO: Check if this can be deleted
        ## The run type is going to be written in the run.log file
        #self.__run_type   = 'basic'

        ## Counters
        ## Number of runs per group
        ## The runs are going to be divided into groups.
        ## For example if we are doing a convergence plot:
        ## One group equals one convergence plot
        ## Usually a group only contains one run
        #self.__no_runs_in_group = False
        ## Count number of runs in a group
        #self.__run_counter      = False
        ## A group counter
        #self.__group_no = 0
        ## Dictionary to be filled with the folders and the status of the
        ## different runs in the  group
        #self.__run_groups = {}
        ## The entries of dmp_folder are the paths where the dmp files
        ## will be sotred. This makes our job easy when we want to make
        ## a plot of grouped runs.
        ## The entries in job_status will be filled with the job_name
        ## If the job has already been done previously, the job status will
        ## be set to 'done'.
        #self.__run_groups[self.__group_no] ={'dmp_folder':[], 'job_status':[]}
        #}}}

        # Check if the instance is set correctly
        self.__check_for_basic_instance_error()
#}}}

# The destructor
# TODO: Get the exit code of the process. If the process failed, write
#       that an error occured instead
#       Could eventually check if a error flag has been set to true
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
            print('\n'*3)
            print(' ' + '~'*69 + '\n'*3)
        elif len(self.__warnings) > 0 and not(error_occured):
            print('\n'*3 + ' ' + '~'*69)
            print("| No WARNINGS detected before instance destruction in"+\
                  " 'bout_runners'. |")
#}}}

# The run function
#{{{run
    def run(self,\
            remove_old = False,\
            post_processing_function = None,\
            post_process_after_every_run = None,\
            **kwargs):
        #{{{docstring
        """
        Makes a run for each of the combination given by the member data.

        Input
        remove_old                  - boolean telling whether old run
                                      files should be deleted or not
        post_processing_function    - a function to be called after
                                      one or several run
        post_process_after_each_run - boolean telling whether
                                      post_processing_function
                                      should be called after each run,
                                      or rather after
        **kwargs                    - parameters to be passed to the
                                      post_processing_function
        """
        #}}}

        # Check for errors in the run function
        self.__error_check_for_run_input(remove_old,\
                                         post_processing_function,\
                                         post_process_after_every_run)

        # TODO: Could possibly delete this then
        # Check for additional errors (if any)
        self.__additional_error_check(**kwargs)

        # Create the run log
        self.__create_run_log()

        # We check that the given combination of nx and ny is
        # possible to perform with the given nproc
        if (self.__nx != None) and (self.__ny != None):
            self.__get_correct_domain_split()

        # Get the combinations of the member functions
        possibilities = self.__get_possibilities()
        combinations = self.__get_combinations(possibilities)

        # Print either 'now running' or 'now submitting'
        self.__print_run_or_submit()

        # TODO: See if you can remove this
        ## Set the run_counter and the number of runs in one group
        #self.__set_run_counter(**kwargs)
        # TODO: Is the run counter superflous as we already have the run
        #       counter?
        # Set the run counter
        self.__run_counter = 0


        # FIXME: THINK THE ANSWER TO OUR PRAYERS ARE HERE - GET THE
        # FOLDER NAMES FROM THE STRINGS :D
        # The run
        for run_no, combination in enumerate(combinations):

            # Get the folder to store the data
            self.__prepare_dmp_folder(combination)

            if remove_old:
                # Remove old data
               self.__remove_data()

            # Check if the run has been performed previously
            do_run = self.__check_if_run_already_performed()
            # FIXME: YOU ARE HERE
            # Do the actual runs
            if do_run:
                # Call the driver for a run
                self.__run_driver(combination, run_no)

                # Add one to the run counter
                self.run_counter +=1

            #{{{ TODO: Check if can delete
            ## Append the current dump folder to the current run_group
            #self.__run_groups[self.__group_no]['dmp_folder'].append(self.__dmp_folder)
            ## Update the run counter
            #self.__run_counter += 1

            ## If we have looped over all the folders in a run
            ## group, we will change the run group (as long as run_no is
            ## lower than all_combinations [remember that run_no starts
            ## from 0 and len() 'starts' from 1])
            #}}}


            # FIXME: Question: would we like to sort, AND NOT do
            #        postprocessing?
            #        Assume yes...because of maybe the run takes long
            #        time or whatever

            # FIXME: YOU ARE HERE
            # If a post processing function is to be given
            # OK, think here, what do we need to pass to the post
            # processing function...the names of all the folder in a
            # group
            # TODO: Let the run function allow kwargs
            #       Let the run function accept a post-processing
            #       function
            if post_processing_after_one_run:
                call the post_processing_function(current_dmp)

            if post_processing_after_group_of_run_finished
                append_the_dmp to dmps
                if...se below:
                    call the post_processing_function(current_dmps)


            # TODO: Check when THE CONSTANT is changing and post a
            #       post-processing commit
            # TODO: Make a flag whether or not to group runs
            if (self.__run_counter % self.__no_runs_in_group == 0)\
                and (run_no < len(all_combinations)-1):

                # Update the group number
                self.__group_no += 1
                # Create a new key in the dictionary for the new
                # run group
                self.__run_groups[self.__group_no] =\
                    {'dmp_folder':[], 'job_status':[]}

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
        output, run_time = self.__single_run( combination )
        # Print info to the log file for the runs
        self.__append_run_log(start, run_no, run_time)
        print('\n')

        # TODO: DELETE THIS
        ## As the jobs are being run in serial, the status of the job
        ## will be done for all jobs
        #self.__run_groups[self.__group_no]['job_status'].append('done')
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
                self.__make()
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
        if self.__MXG != None and type(self.__MXG) != int:
            message  = "MXG is of wrong type\n"+\
                       "MXG must be given as an int"
            self.__errors.append("TypeError")
            raise TypeError(message)
        if self.__MYG != None and type(self.__MYG) != int:
            message  = "MYG is of wrong type\n"+\
                       "MYG must be given as an int"
            self.__errors.append("TypeError")
            raise TypeError(message)
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

        #{{{Check if nx, ny, nz and nout are int or list of int
        # Check if the following is an integer, or an iterable
        # containing only integers
        check_if_int = [\
            (self.__nx        , 'nx')        ,\
            (self.__ny        , 'ny')        ,\
            (self.__nz        , 'nz')        ,\
            (self.__nout      , 'nout')      ,\
            ]

        self.__check_for_correct_type(var = check_if_int,\
                                    the_type = int)
        #}}}

        #{{{Check if timestep is Number or list of Number
        # Check if the following is a number
        check_if_number = [\
            (self.__timestep  , 'timestep')\
            ]

        self.__check_for_correct_type(var = check_if_number,\
                                    the_type = Number)
        #}}}

        #{{{Check if solver, grid_file and methos is str or list of str
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
            ]

        self.__check_for_correct_type(var = check_if_string,\
                                    the_type = str)
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

        #{{{Check if sort_by is set correctly
        if self.__sort_by != None:
            if type(self.__sort_by) != str:
                self.__errors.append("TypeError")
                raise TypeError ("sort_by must be set as a string when set")

        possible_method = [\
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

        # Append the additionals
        # If additional is set
        if self.__additional != None:
            for additional in self.__additional:
                # The additional now contains a tuple of three elements
                # We would like to extract the name of them, and append
                # it to the possibilities list
                possible_method.append(additional[0])

        # Make a list of the variables
        the_vars = [\
            (self.__sort_by, 'sort_by')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] != None:
                self.__check_if_set_correctly(var           = var,\
                                              possibilities = possible_method)
        #}}}

        #{{{Check len of nx, ny, nz and that consistent with grid_file
        # Check if nx, ny, nz and gridfile is set at the same time
        if (self.__nx != None and self.__grid_file != None) or\
           (self.__ny != None and self.__grid_file != None) or\
           (self.__nz != None and self.__grid_file != None):
            # Read the variable from the file
            for grid_file in self.__grid_file:
                # Open (and automatically close) the grid files
                f = DataFile(grid_file)
                # Search for nx, ny and nz in the grid file
                domain_types = ["nx", "ny", "nz"]
                for domain_type in domain_types:
                    grid_variable = f.read(domain_type)
                    # If the variable is found
                    if grid_variable != None:
                        self.__errors.append("TypeError")
                        message  = domain_type + " was specified both in the "
                        message += "driver and in the grid file.\n"
                        message += "Please remove " + domain_type
                        message += " from the driver if you would "
                        message += "like to run with a grid file."
                        raise TypeError(message)

        # If grid files are set, use the nx, ny and nz values in the
        # member data if applicable
        if self.__grid_file != None:
            # Make a dict of appendable lists
            spatial_domain = {'nx':[], 'ny':[], 'nz':[]}
            for grid_file in self.__grid_file:
                # Open (and automatically close) the grid files
                f = DataFile(grid_file)
                # Search for nx, ny and nz in the grid file
                domain_types = ["nx", "ny", "nz"]
                for domain_type in domain_types:
                    grid_variable = f.read(domain_type)
                    # If the variable is found
                    if grid_variable != None:
                        spatial_domain[domain_type].append(grid_variable)
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

        # Check that nx, ny and nz are set correctly
        if self.__nx != None and self.__ny != None:
            self.__check_if_same_len((self.__nx, 'nx'), (self.__ny, 'ny'))
        if self.__nx != None and self.__nz != None:
            self.__check_if_same_len((self.__nx, 'nx'), (self.__nz, 'nz'))
        if self.__ny != None and self.__nz != None:
            self.__check_if_same_len((self.__ny, 'ny'), (self.__nz, 'nz'))
        #}}}

        #{{{ Check that timestep and nout have the same len
        if self.__timestep != None and self.__nout != None:
            self.__check_if_same_len((self.__timestep, 'timestep'),\
                                   (self.__nout, 'nout'))
        #}}}

        #{{{Check that additional is on the correct form
        # additional should be on the form
        # additional = [(section1, name1, [value1-1, value1-2, ...]),\
        #               (section2, name2, [value2-1, value2-2, ...]),\
        #               ...]
        # We will now check that
        if self.__additional != None:
            # Set a success variable that will fail if anything goes
            # wrong
            success = True
            # Check if self.__addition is iterable
            if hasattr(self.__additional, "__iter__"):
                # Check if self.__additional is a string
                if type(self.__additional) != str and\
                   type(self.__additional) != dict:
                    # If additional is given as an iterable
                    if hasattr(self.__additional[0], "__iter__" ):
                        # Do the same check as above for all the
                        # elements
                        for elem in self.__additional:
                            # Check if self.__addition is iterable
                            if hasattr(elem, "__iter__"):
                                # Check if elem is a string
                                if type(elem) != str:
                                    if type(elem[0]) == str:
                                        # Check that the second element
                                        # (the name) is a string
                                        if type(elem[1]) != str:
                                            success = False
                                        # If more than three elements
                                        # are given
                                        if len(elem) != 3:
                                            success = False
                                    # elem[0] is not a string
                                    else:
                                        success = False
                                # elem is a string
                                else:
                                    success = False
                            # elem is not iterable
                            else:
                                success = False
                    # self.__additional[0] is not a string, and not iterable
                    else:
                        success = False
                # self.__additional is a string or a dict
                else:
                    success = False
            # self.__additional is not iterable
            else:
                success = False
            if not(success):
                message  = "self.__additional is on the wrong form.\n"
                message += "self.__additional should be on the form\n"
                message += "self.__additional=\ \n"
                message +=\
                        "     [(section1, name1, [value1-1, value1-2,...]),\ \n"
                message +=\
                        "      (section2, name2, [value2-1, value2-2,...]),\ \n"
                message +=\
                        "       ...])\n"
                self.__errors.append("TypeError")
                raise TypeError(message)
        #}}}

        #{{{Check if cpy_source is boolean
        check_if_bool = [\
            (self.__cpy_source, 'cpy_source'),\
            (self.__allow_size_modification, 'allow_size_modification')\
            ]

        self.__check_for_correct_type(var = check_if_bool,\
                                    the_type = bool)
        #}}}
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

        #{{{Check if remove_old and restart is set on the same time
        if remove_old == True and self.__restart != None:
            self.__errors.append("RuntimeError")
            raise RuntimeError("You should not remove old data if you"\
                               " want a restart run")
        #}}}

        #{{{Check that the post_processing_function is a fuction
        if (post_processing_function != None) and\
           (type(post_processing_function) != function):
            self.__errors.append("RuntimeError")
            message = "post_process_after_every_run must be a"+\
                      " function"
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

# TODO: Could probably delete this as we are now doing the error check
#       in the constructor
#{{{__additional_error_check
    def __additional_error_check(self, **kwargs):
        """Virtual function. Will in child classes check for additional
        errors"""
        return
#}}}

#{{{__create_run_log
    def __create_run_log(self):
        """Makes a run_log file if it doesn't exists"""

        # Checks if run_log exists
        self.__run_log = self.__directory + "/run_log.txt"
        if os.path.isfile(self.__run_log) == False:
            # Create a file to be appended for each run
            f = open(self.__run_log , "w")
            # The header
            header = ['start_time', 'run_type', 'run_no', 'dump_folder', 'run_time_H:M:S']
            header = '    '.join(header)
            f.write('#' + header + '\n')
            f.close()

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
            self.__MXG = self.__find_variable_in_BOUT_inp("MXG")

            # Check that one unique variable was found
            if len(self.__MXG) > 1:
                self.__errors.append("RuntimeError")
                message =  "Several matches was found when searching for "
                message += "MXG.\n"
                message += "Please check " + self.__directory + "/BOUT.inp"
                raise RuntimeError(message)

            # Pick the only element
            self.__MXG = int(self.__MXG[0])

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
                    MX = self.__nx[size_nr] - 2*self.__MXG
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
                self.__warnings.append(message)
#}}}

#{{{__get_possibilities
    def __get_possibilities(self):
        """ Returns the list of the possibilities. In get_combinations
        the elements of this list is going to be put together to a list
        of strings which will be used when making a run."""

        # Set the combination of nx, ny and nz (if it is not already
        # given by the gridfile)
        # Appendable list
        spatial_grid_possibilities = []
        if (self.__grid_file == None):
            # Appendable lists
            nx_str = []
            ny_str = []
            nz_str = []
            # Append the different dimension to the list of strings
            if self.__nx != None:
                for nx in self.__nx:
                    nx_str.append('mesh:nx=' + str(nx))
            if self.__ny != None:
                for ny in self.__ny:
                    ny_str.append('mesh:ny=' + str(ny))
            if self.__nz != None:
                for nz in self.__nz:
                    nz_str.append('mesh:nz=' + str(nz))
            # Combine the strings to one string
            # Find the largest length
            max_len = np.max([len(nx_str), len(ny_str), len(nz_str)])
            # Make the strings the same length
            if len(nx_str) < max_len:
                nx_str.append('')
            if len(ny_str) < max_len:
                ny_str.append('')
            if len(nz_str) < max_len:
                nz_str.append('')
            # Append the spatial grid possibilities as a string
            for number in range(max_len):
                # Make a list
                current_grid = [nx_str[number],\
                                ny_str[number],\
                                nz_str[number]]
                # Join the strings in the list and append
                spatial_grid_possibilities.append(' '.join(current_grid))

        # Set the combination of timestep and nout if set
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

        # List of the possibilities of the different variables
        list_of_possibilities = [spatial_grid_possibilities,\
                                 temporal_grid_possibilities]

        # Put MXG and MYG into a list if they are not set to None
        # This makes the memberdata iterable, and useable in
        # generate_possibilities
        if self.__MXG != None:
            self.__MXG = [self.__MXG]
        if self.__MYG != None:
            self.__MYG = [self.__MYG]

        # List of tuple of varibles to generate possibilities from
        tuple_of_variables = [\
            (self.__solver,     "solver", "type"),\
            (self.__grid_file,  "",       "grid"),\
            (self.__ddx_first,  "ddx",    "first"),\
            (self.__ddx_second, "ddx",    "second"),\
            (self.__ddx_upwind, "ddx",    "upwind"),\
            (self.__ddx_flux,   "ddx",    "flux"),\
            (self.__ddy_first,  "ddy",    "first"),\
            (self.__ddy_second, "ddy",    "second"),\
            (self.__ddy_upwind, "ddy",    "upwind"),\
            (self.__ddy_flux,   "ddy",    "flux"),\
            (self.__ddz_first,  "ddz",    "first"),\
            (self.__ddz_second, "ddz",    "second"),\
            (self.__ddz_upwind, "ddz",    "upwind"),\
            (self.__ddz_flux,   "ddz",    "flux"),\
            (self.__MXG,        "",       "MXG"),\
            (self.__MYG,        "",       "MYG"),\
            ]

        # Append the additional option to tuple of variables
        for additional in self.__additional:
            # If the last element of additional is not iterable we need
            # put them into a list to make them iterable (in order to
            # use them in generate_possibilities)
            if not(hasattr(additional[2], "__iter__")):
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

        # Append the possibilities to the list of possibilities
        for var in tuple_of_variables:
            list_of_possibilities.append(\
                    self.__generate_possibilities(var[0], var[1], var[2])\
                    )

        # Return the list_of possibilities
        return list_of_possibilities
#}}}

# TODO: SHOULD ADD THE OPTION: FASTEST VARYING
#{{{__get_combinations
    def __get_combinations(self, input_list):
        """ Takes a list with lists as element as an input.
        Returns a list of all combinations between the elements of the
        lists of the input list """

        # Remove empty elements in input_list in order for
        # itertools.product to work
        input_list = [elem for elem in input_list if elem != []]

        # TODO: Rewrite the comments
        # The list of combinations is a list with tuples as elements
        # We would like to combind the element in these tuples to one
        # string
        # The last element in the list will be the fastest varying
        # element
        print(input_list)
        import pdb
        pdb.set_trace()
        # TODO: GET THE STRING FROM THIS FASTEST VARYING?

        # Sort the input list (so that the fastest varying variable can
        # be chosen)

        # If we would like to sort
        if self.__sort_by != None:
            # Swap the list corresponding to the sort_by statement so
            # that that list will be the last. The itertools.product
            # will then make that list the fastest varying in the list
            input_list = self.__get_swapped_input_list()

        all_combinations_as_tuple = list(itertools.product(*input_list))

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
        # Obtain file and folder names
        folder_name = self.__get_folder_name(combination)
        # Create folder if it doesn't exists
        self.__dmp_folder = self.__directory + "/" + folder_name
        create_folder(self.__dmp_folder)
        # Copy the input file into this folder
        command = 'cp ' + self.__directory + '/BOUT.inp ' + self.__dmp_folder
        shell(command)

        # Copy the source files if cpy_source is True
        self.__cpy_source = True
        if self.__cpy_source:
            # This will copy all C++ files to the dmp_folder
            cpp_extension= ['.cc', '.cpp', '.cxx', '.C', '.c++',\
                            '.h',  '.hpp', '.hxx', '.h++']
            # Copy for all files in the extension
            for extension in cpp_extension:
                file_names = glob.glob('*' + extension)
                for a_file in file_names:
                    command = 'cp ' + a_file + ' ' + self.__dmp_folder
                    shell(command)
#}}}

#{{{__remove_data
    def __remove_data(self):
        """Removes *.nc and *.log files from the dump directory"""

        print("Removing old data")
        # Make the command
        command = "rm -f ./" + self.__dmp_folder +\
                  "/*.nc ./" + self.__dmp_folder +\
                  "/*.log"
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
            print('To overwrite old files, run with self.__run(remove_old=True)\n')
            return False
        else:
            return True
#}}}

#{{{__post_run
    def __post_run(self, **kwarg):
        """In basic_runner this is a virtual function"""
        return
#}}}
#}}}

#{{{Function called by __set_program_name
#{{{__make
    def __make(self):
        """Makes the .cxx program, saves the make.log and make.err"""
        print("Making the .cxx program")
        command = "make > make.log 2> make.err"
        shell(command)
        # Check if any errors occured
        if os.stat("make.err").st_size != 0:
            self.__errors.append("RuntimeError")
            raise RuntimeError("Error encountered during make, see 'make.err'.")
#}}}
#}}}

#{{{ Functions called by the basic_error_checker
#{{{__check_for_correct_type
    def __check_for_correct_type(self,\
                                 var      = None,\
                                 the_type = None):
        """Checks if a varible has the correct type

        Input:
        var      - a tuple consisting of
                   var[0] - the variable (a data member)
                   var[1] - the name of the varible given as a string
        the_type - the data type"""

        # Set a variable which is False if the test fails
        success = True
        for cur_var in var:
            # There is an option that the variable could be set to None,
            # and that the default value from BOUT.inp will be used
            if cur_var[0] != None:
                # Check for the correct type
                if isinstance(cur_var[0], the_type) == False:
                    # Check if it is an iterable
                    if hasattr(cur_var[0], "__iter__") and\
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
                               cur_var[1] + " must be " + the_type.__name__  +\
                               " or an iterable with " + the_type.__name__ +\
                               " as elements."
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
        form solver/method/additional/ghost_nout_timestep/mesh"""

        # Combination is one of the combination of the data members
        # which is used as the command line arguments in the run
        combination = combination.split()

        # Make lists for the folder-type, so that we can append the
        # elements in the combination folders if it is found
        solver              = []
        method              = []
        additional          = []
        ghost_nout_timestep = []
        mesh                = []

        # We will loop over the names describing the methods used
        # Possible directional derivatives
        dir_derivatives = ['ddx', 'ddy', 'ddz']

        # Check trough all the elements of combination
        for elem in combination:
            # If solver is in the element
            if 'solver' in elem:
                # Remove ':type', and append it to the
                # solver folder
                cur_solver = elem.replace(':type','')
                cur_solver = cur_solver.replace('solver=','')
                solver.append(cur_solver)
            # If MXG, MYG, nout or timestep is in the element
            elif ('MXG' in elem) or\
               ('MYG' in elem) or\
               ('nout' in elem) or\
               ('timestep' in elem):
                # Remove '=', and append it to the
                # ghost_nout_timestep folder
                ghost_nout_timestep.append(elem.replace('=','_'))
            # If nx, ny or nz is in the combination
            elif 'mesh' in elem:
                # Remove 'mesh:', and append it to the mesh folder
                cur_mesh = elem.replace('mesh:','')
                cur_mesh = cur_mesh.replace('=','_')
                mesh.append(cur_mesh)
            # If the element is none of the above
            else:
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
        method.sort()
        additional.sort()
        ghost_nout_timestep.sort()
        mesh.sort()

        # Combine the elements in the various folders
        method              = ['_'.join(method)]
        additional          = ['_'.join(additional)]
        ghost_nout_timestep = ['_'.join(ghost_nout_timestep)]
        mesh                = ['_'.join(mesh)]

        # Put all the folders into the combination_folder
        # We access the zeroth element as the folders are given as a
        # list
        combination_folder = [solver[0],\
                              method[0],\
                              additional[0],\
                              ghost_nout_timestep[0],\
                              mesh[0]]

        # Make the combination folder as a string
        combination_folder = '/'.join(combination_folder)

        return combination_folder
#}}}
#}}}

#{{{Function called by __run_driver
#{{{__single_run
    def __single_run(self, combination=''):
        """Makes a single MPIRUN of the program"""

        # Get the command to be used
        command = self.__get_command_to_run( combination )

        import pdb
        pdb.set_trace()
        # Time how long the time took
        tic = timeit.default_timer()

        # Launch the command
        status, out = launch(command,\
                             runcmd = self.__MPIRUN,\
                             nproc = self.__nproc,\
                             pipe = True,\
                             verbose = True)

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
        line = [start_time, self.__run_type, run_no, dmp_line, run_time]
        # Opens for appending
        f = open(self.__run_log , "a")
        f.write('    '.join(str(element) for element in line) + "\n")
        f.close()
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

        # Initialize the length of a group (see below for explanaition)
        self.__len_group = None
        # Find what list in the list which contains what we would sort
        # by
        # If we would like to sort by the spatial domain
        if self.__sort_by == 'spatial_domain':
            # nx, ny and nz are all under the section 'mesh'
            find_in_list = 'mesh'
            # Text to be printed if find_in_list is not found
            text_if_not_found = 'neither nx, ny or nz was found'
        # If we would like to sort by the temporal domain
        elif self.__sort_by == 'temporal_domain':
            # If we are sorting by the temporal domain, we can either
            # search for timestep or nout
            if self.__timestep != None:
                find_in_list = 'timestep'
            elif self.__nout != None:
                find_in_list = 'nout'
            else:
                message = "Could not sort by 'temporal_domain' as"+\
                          " neither 'timestep' nor 'nout' is given")
                raise RuntimeError(message)
        # If we would like to sort by the method
        elif (self.__sort_by == 'ddx_first' or\
             (self.__sort_by == 'ddx_second') or\
             (self.__sort_by == 'ddx_upwind') or\
             (self.__sort_by == 'ddx_flux') or\
             (self.__sort_by == 'ddy_first') or\
             (self.__sort_by == 'ddy_second') or\
             (self.__sort_by == 'ddy_upwind') or\
             (self.__sort_by == 'ddy_flux') or\
             (self.__sort_by == 'ddz_first') or\
             (self.__sort_by == 'ddz_second') or\
             (self.__sort_by == 'ddz_upwind') or\
             (self.__sort_by == 'ddz_flux'):
            find_in_list = self.__sort_by.replace('_',':')
            text_if_not_found = self.__sort_by + " was not found"
        # One is either sorting by solver
        elif self.__sort_by == 'solver':
            find_in_list = self.__sort_by
            text_if_not_found = self.__sort_by + "was not found"
        # One is solving by additional
        else:
            text_if_not_found = self.__sort_by + " was not given as"+\
                                " an 'additional' option"


        # Loop over the lists in the input list
        # Make a flag to break the outermost loop if find_in_list is
        # found
        break_outer = False
        for elem_nr, elem in enumerate(input_list):
            # Each of the elements in this list is a string
            for string in elem:
                # Check if fins_in_list is in the string
                if find_in_list in string:
                    # If there is a match, store the element number
                    swap_from_index = elem_nr
                    # Check the length of the element (as this is the
                    # number of times the run is repeated, only changing
                    # the values of sort_by [defining a group])
                    self.__len_group = len(elem)
                    # Break the loop to save time
                    break_outer = True
                    break
            # Break the outer loop if find_in_list_is_found
            if break_outer:
                break

        # If there was no match, throw an error
        if self.__len_group == None:
            message  = "Could not sort by " + self.__sort_by + " as "
            message += text_if_not_found
            raise RuntimeError(message)

        # As it is the last index which changes the fastest, we swap the
        # element where the find_in_list was found with the last element
        input_list[swap_from_index], input_list[-1] =\
                input_list[-1], input_list[swap_from_index]

        return input_list
#}}}
#}}}

#{{{Function called by __single_run
#{{{__get_command_to_run
    def __get_command_to_run(self, combination):
        """ Returns a string of the command which will run the BOUT++
        program"""

        # Creating the arguments
        arg = " -d " + self.__dmp_folder + combination

        # If the run is restarted with initial values from the last run
        if self.__restart != False:
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

# TODO: Check if this can be deleted
#{{{__set_run_counter
    def __set_run_counter(self, **kwargs):
        """Sets self.__run_counter and self.__no_runs_in_group"""

        self.__run_counter = 0
        # Normally there is only one job in a group
        self.__no_runs_in_group = 1

        # Due to plotting it is convenient to put runs belonging to the
        # same convergence plot into one group
        if ('convergence_type' in list(kwargs.keys())):
            if kwargs['convergence_type'] == 'spatial':
                # How many runs must we make before we can make a
                # convergence plot
                keys = list(self.__n_points.keys())
                self.__no_runs_in_group = len(self.__n_points[keys[0]])
            elif kwargs['convergence_type'] == 'temporal':
                # How many runs must we make before we can make a
                # convergence plot
                self.__no_runs_in_group = len(self.__timestep)
#}}}
#}}}



#{{{basic_qsub_runner
class basic_qsub_runner(basic_runner):
#{{{docstring
    """Class for running BOUT++.
    Works like the basic_runner, but submits the jobs to a torque queue
    with qsub.

    The link below gives a nice introduction to the qsub system
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
        # 'super' from qsub_run_with_plots, which inherits from both
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
        self.__run_type   = 'basic_qsub'
        # A string which will be used to write a self.__deleting python
        # script
        self.__python_tmp = ''
        # The jobid returned from the qsub
        self.__qsub_id = None
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
        self.__qsub_error_check(**kwargs)
#}}}

#{{{single_submit
    def single_submit(self, combination, run_no):
        """Single qsub submission"""
        # Get the script (as a string) which is going to be
        # submitted
        job_name, job_string =\
            self.__get_job_string(run_no, combination)

        # The submission
        self.__qsub_id = self.__submit_to_qsub(job_string)
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
        # qsub in
        create_folder(self.__directory + '/qsub_output')

        # Get the start_time
        start_time = self.__get_start_time()

        # The name of the file
        python_name = 'clean_up_'+start_time+'.py'

        # Creating the job string
        job_name = 'clean_up_' + self.__run_type + '_'+ start_time

        # Get the core of the job_string (note that we only need to use
        # one node and one processor for this)
        job_string = self.__create_qsub_core_string(\
            job_name, '1', '1', self.__walltime,\
            folder = self.__directory + '/qsub_output/')
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
        self.__submit_to_qsub(job_string, dependent_job = self.__qsub_id)
#}}}
#}}}

# Auxiliary functions
#{{{
#{{{qsub_error_check
    def qsub_error_check(self, **kwargs):
        """Checks for specific qsub errors"""
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
                                     " qsub runners must be strings")
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
        qsub."""

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
        job_string = self.__create_qsub_core_string(\
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

#{{{create_qsub_core_string
    def create_qsub_core_string(\
        self, job_name, nodes, ppn, walltime, folder=''):
        """Creates the core of a qsub script as a string"""

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

#{{{submit_to_qsub
    def submit_to_qsub(self, job_string, dependent_job=None):
        """Saves the job_string as a shell script, submits it and
        deletes it. Returns the output from qsub as a string"""
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



# TODO: Delete this
#{{{class run_with_plots
class run_with_plots(basic_runner):
#{{{docstring
    """Class running BOUT++ in the same way as the basic_runner, with
    the additional feature that it calls one of the plotters in
    'bout_plotters'.

    For further details, see the documentation of basic_runner.
    """
#}}}

# The constructor
#{{{__init__
    def __init__(self,\
                 plot_type  = False,\
                 extension  = 'png',\
                 solver     = False,\
                 nproc      = 1,\
                 methods    = False,\
                 n_points   = False,\
                 directory  = 'data',\
                 nout       = False,\
                 timestep   = False,\
                 MXG        = False,\
                 MYG        = False,\
                 additional = False,\
                 restart    = False,\
                 **kwargs):
        """Specify either how many time indices you want to plot in one
        plot, or give a list of time indices. If both are given, then the
        list of plot_times has higher precedence."""

        # Note that the constructor accepts additional keyword
        # arguments. This is because the constructor can be called with
        # 'super' from qsub_run_with_plots, which inherits from both
        # basic_qsub_runner and run_with_plots (which takes different
        # arguments as input)

        # Call the constructor of the superclass
        super(run_with_plots,  self.__init__(solver     = solver,\
                                             nproc      = nproc,\
                                             methods    = methods,\
                                             n_points   = n_points,\
                                             directory  = directory,\
                                             nout       = nout,\
                                             timestep   = timestep,\
                                             MXG        = MXG,\
                                             MYG        = MYG,\
                                             additional = additional,\
                                             restart    = restart,\
                                             **kwargs))

        if plot_type == False:
            self.__errors.append("TypeError")
            raise TypeError ("Keyword argument 'plot_type' must be given"+\
                             " when running run_with_plots")

        self.__plot_type           = plot_type
        self.__file_extension      = extension
        self.__run_type            = 'plot_' + self.__plot_type

        # Check if a DISPLAY is set
        try:
            os.environ['DISPLAY']
        except KeyError:
            message =  "No display is set! Changing the backend to 'Agg' in"+\
                       " order to plot."
            self.__warnings.append(message)
            warning_printer(message)
            import matplotlib.pyplot as plt
            plt.switch_backend('Agg')
#}}}

# Functions called directly by the main function
#{{{
#{{{additional_error_check
    def additional_error_check(self, **kwargs):
        """Checks for errors related to the relevant plotter to the current
        class"""

        # Since the error checkers should be called from both
        # bout_runners and bout_plotters, the error_checkers have been
        # put in the class check_for_plotters_errors defined in
        # common_bout_functions
        # The error checker is called by the constructor, so all we have
        # to do is to create an instance of the class
        plotter_error_checker =\
           check_for_plotters_errors(self.__plot_type, n_points=self.__n_points,
                                     timestep=self.__timestep, **kwargs)
#}}}

#{{{post_run
    def post_run(self, **kwargs):
        """Calls self.__plotter_chooser"""
        self.__plotter_chooser(**kwargs)
#}}}
#}}}

# Plotter specific
#{{{
#solution_plot specific
#{{{
#{{{solution_plotter
    def solution_plotter(self,\
            show_plots = False,\
            collect_x_ghost_points = False, collect_y_ghost_points = False,\
            **kwargs):
        """Calls the correct plotter from bout_plotters"""

        # Creates an instance of solution_plotter
        make_my_sol_plot = solution_plotter(\
            run_groups             = self.__run_groups,\
            directory              = self.__directory,\
            file_extension         = self.__file_extension,\
            show_plots             = show_plots,\
            collect_x_ghost_points = collect_x_ghost_points,\
            collect_y_ghost_points = collect_y_ghost_points,\
            variables              = kwargs['variables'],\
            plot_direction         = kwargs['plot_direction'],\
            plot_times             = kwargs['plot_times'],\
            number_of_overplots    = kwargs['number_of_overplots'])

        # Run the plotter
        make_my_sol_plot.collect_and_plot()
#}}}
#}}}

#solution_and_error_plotter specific
#{{{
#{{{solution_and_error_plotter
    def solution_and_error_plotter(self,\
            show_plots = False,\
            collect_x_ghost_points = False, collect_y_ghost_points = False,\
            **kwargs):
        """Calls the correct plotter from bout_plotters"""

        # Creates an instance of solution_and_error_plotter
        make_my_sol_err_plot = solution_and_error_plotter(\
            run_groups             = self.__run_groups,\
            directory              = self.__directory,\
            file_extension         = self.__file_extension,\
            show_plots             = show_plots,\
            collect_x_ghost_points = collect_x_ghost_points,\
            collect_y_ghost_points = collect_y_ghost_points,\
            variables              = kwargs['variables'],\
            plot_direction         = kwargs['plot_direction'],\
            plot_times             = kwargs['plot_times'],\
            number_of_overplots    = kwargs['number_of_overplots'])

        # Run the plotter
        make_my_sol_err_plot.collect_and_plot()
#}}}

#{{{get_plot_times_and_number_of_overplots
    def get_plot_times_and_number_of_overplots(self, **kwargs):
        """Returns plot_times and number_of_overplots"""

        if self.__plot_type == 'solution_and_error_plot' or\
           self.__plot_type == 'solution_plot':
            kwarg_keys = list(kwargs.keys())

            if ('plot_times' in kwarg_keys) == False:
                plot_times = False
            else:
                plot_times = kwargs['plot_times']
            if ('number_of_overplots' in kwarg_keys) == False:
                number_of_overplots = False
            else:
                number_of_overplots = kwargs['number_of_overplots']

        return plot_times, number_of_overplots
#}}}
#}}}

#convergence_plot specific
#{{{
#{{{convergence_plotter
    def convergence_plotter(self,\
            show_plots = False,\
            collect_x_ghost_points = False, collect_y_ghost_points = False,\
            **kwargs):
        """Calls the convergence plotter from bout_plotters"""
        # Creates an instance of convergence_plotter
        make_my_convergence_plots = convergence_plotter(\
           run_groups             =  self.__run_groups,\
           directory              =  self.__directory,\
           file_extension         =  self.__file_extension,\
           show_plots             =  show_plots,\
           collect_x_ghost_points =  collect_x_ghost_points,\
           collect_y_ghost_points =  collect_y_ghost_points,\
           variables              =  kwargs['variables'],\
           convergence_type       =  kwargs['convergence_type'])

        # Run the plotter
        make_my_convergence_plots.collect_and_plot()
#}}}

#{{{set_nx_range
    def set_nx_range(self, grid_keys, inner_points, ranges):
        """Append ranges (a list filled with the number of grid points) with
        the grid points in nx (given from the list 'inner ranges')"""
        # We must find the guard cells in x to find the
        # correct nx
        # We search first if MXG is set
        if self.__MXG:
            MXG = self.__MXG
        else:
            # We must find it in the input file
            MXG = find_variable_in_BOUT_inp(self.__directory,\
                                            'MXG')
            # If MXG was not found
            if type(MXG) == str:
                # MXG is set to the default value
                MXG = 2
        # Write the range to a list of strings
        # Set the x_range
        x_range = ['mesh:nx=' + str(nr+2*MXG) for nr in \
                   inner_points]
        ranges.append(x_range)
        return ranges
#}}}

#{{{set_ny_range
    def set_ny_range(self, grid_keys, inner_points, ranges):
        """Append ranges (a list filled with the number of grid points) with
        the grid points in ny (given from the list 'inner_points')"""
        # ny is the number of inner_points
        # Write the range to a list of strings
        # Set the y_range
        y_range = ['mesh:ny=' + str(nr) for nr in inner_points]
        ranges.append(y_range)
        return ranges
#}}}

#{{{set_nz_range
    def set_MZ_range(self, grid_keys, inner_points, ranges):
        """Append ranges (a list filled with the number of grid points) with
        the grid points in MZ (given from the list 'inner_points')"""
        # MZ is the number of inner_points + 1
        # Write the range to a list of strings
        # Set the z_range
        z_range = ['MZ=' + str(nr+1) for nr in inner_points]
        ranges.append(z_range)
        return ranges
#}}}
#}}}
#}}}

# Auxiliary functions
#{{{
#{{{plotter_chooser
    def plotter_chooser(self, **kwargs):
        """Calls the correct plotter from bout_plotters"""

        # Set option
        kwarg_keys = list(kwargs.keys())
        if 'show_plots' in kwarg_keys:
            show_plots = kwargs['show_plots']
        else:
            show_plots = False
        if 'collect_x_ghost_points' in kwarg_keys:
            collect_x_ghost_points = kwargs['collect_x_ghost_points']
        else:
            collect_x_ghost_points = False
        if 'collect_y_ghost_points' in kwarg_keys:
            collect_y_ghost_points = kwargs['collect_y_ghost_points']
        else:
            collect_y_ghost_points = False

        if self.__plot_type == 'solution_and_error_plot' or\
           self.__plot_type == 'solution_plot':
            # Get the plot_times or the number_of_overplots (one may
            # have to be set to False if not given by the user input)
            plot_times, number_of_overplots =\
                self.__get_plot_times_and_number_of_overplots(**kwargs)

            if self.__plot_type == 'solution_plot':
                # Call the solution plotter
                self.__solution_plotter(\
                    plot_times = plot_times,\
                    number_of_overplots = number_of_overplots,\
                    show_plots = show_plots,\
                    collect_x_ghost_points = collect_x_ghost_points,\
                    collect_y_ghost_points = collect_y_ghost_points,\
                    variables = kwargs['variables'],\
                    plot_direction = kwargs['plot_direction'])
            elif self.__plot_type == 'solution_and_error_plot':
                # Call the solution and error plotter
                self.__solution_and_error_plotter(\
                    plot_times = plot_times,\
                    number_of_overplots = number_of_overplots,\
                    show_plots = show_plots,\
                    collect_x_ghost_points = collect_x_ghost_points,\
                    collect_y_ghost_points = collect_y_ghost_points,\
                    variables = kwargs['variables'],\
                    plot_direction = kwargs['plot_direction'])
        elif self.__plot_type == 'convergence_plot':
            # Call the convergence_plotter
            self.__convergence_plotter(\
                show_plots = show_plots,\
                collect_x_ghost_points = collect_x_ghost_points,\
                collect_y_ghost_points = collect_y_ghost_points,\
                variables = kwargs['variables'],\
                convergence_type = kwargs['convergence_type'])
        else:
            self.__errors.append("TypeError")
            raise TypeError ("The given 'plot_type' '" + str(self.__plot_type) +\
                             "' is invalid. See run_with_plots"+\
                             " documentation for valid possibilities.")
#}}}
#}}}
#}}}



# TODO: Delete this
#{{{qsub_run_with_plots
# Note that the basic_qsub_runner runner is inherited first, since we
# want functions in this class to have precedence over the ones run_with_plots
class qsub_run_with_plots(basic_qsub_runner, run_with_plots):
#{{{docstring
    """Class running BOUT++ in the same way as the basic_qsub_runner, with
    the additional feature that it calls one of the plotters in
    'bout_plotters'.

    For further details, see the documentation of basic_qsub_runner and
    run_with_plots.
    """
#}}}

# The constructor
#{{{__init__
    def __init__(self,\
                 plot_type  = False,\
                 extension  = 'png',\
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
                 MXG        = False,\
                 MYG        = False,\
                 additional = False,\
                 restart    = False):
        """Calls the constructors of the super classes."""
        # The cleanest way of doing this is explained in the following link
        # http://stackoverflow.com/questions/13124961/how-to-pass-arguments-efficiently-kwargs-in-python
        # Call the constructors of the super classes
        super(qsub_run_with_plots, self.__init__(plot_type  = plot_type,\
                                                  extension  = extension,\
                                                  nodes      = nodes,\
                                                  ppn        = ppn,\
                                                  walltime   = walltime,\
                                                  mail       = mail,\
                                                  queue      = queue,\
                                                  solver     = solver,\
                                                  nproc      = nproc,\
                                                  methods    = methods,\
                                                  n_points   = n_points,\
                                                  directory  = directory,\
                                                  nout       = nout,\
                                                  timestep   = timestep,\
                                                  MYG        = MYG,\
                                                  MXG        = MXG,\
                                                  additional = additional,\
                                                  restart    = restart))

        self.__run_type = 'qsub_plot_' + self.__plot_type
##}}}

# Functions called directly by the main function
#{{{
#{{{additional_error_check
    def additional_error_check(self, **kwargs):
        """Calls all the error checkers"""

        # Call the error checker for plotter errors from
        # common_bout_functions
        plotter_error_checker =\
           check_for_plotters_errors(self.__plot_type, n_points=self.__n_points,
                                     timestep=self.__timestep, **kwargs)

        # Call the qsub error checker
        self.__qsub_error_check()
#}}}

#{{{post_run
# This should not be a post_run at all, but a common for all plotters
    def post_run(self, **kwargs):
        """Submits a job which calls the plotter from bout_plotters"""

        # Make folder to move the error and log files
        create_folder(self.__directory + '/qsub_output')

        # Get the start_time
        start_time = self.__get_start_time()

        # The name of the file
        python_name = 'tmp'+start_time+'.py'

        # Creating the job string
        job_name = 'post_' + self.__run_type + '_'+ start_time

        # The wall-time of the post-processing should be longer than the
        # wall-time of the runs, as the post-processing will wait for
        # the runs to finish
        # Since plotting should not take to long, we add one hour to the
        # wall time of the run
        # Convert the string to a list of strings
        walltime = self.__walltime.split(':')
        # Convert them to a string
        walltime = [int(element) for element in walltime]
        # Add an hour to the 'hours'
        walltime[0] += 1
        # Convert walltime back to strings
        walltime = [str(element) for element in walltime]
        # Check that the lenght is two (format needs to be HH:MM:SS)
        for nr in range(len(walltime)):
            if len(walltime[nr]) < 2:
                walltime[nr] = '0' + walltime[nr]
        # Join the list of strings
        walltime = ':'.join(walltime)

        # Get the core of the job_string
        job_string = self.__create_qsub_core_string(\
            job_name, nodes='1', ppn='1',\
            walltime=walltime, folder=self.__directory + '/qsub_output/')

        # We will write a python script which calls the
        # relevant bout_plotter

        # First line of the script
        self.__python_tmp = 'import os\n'

        # Append self.__python_tmp with the right plotter
        self.__plotter_chooser(**kwargs)

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
        print('\nSubmitting a job which calls ' + self.__plot_type)
        self.__submit_to_qsub(job_string)
#}}}
#}}}


# Plotter specific
#{{{
#solution_plot specific
#{{{
#{{{solution_plotter
    def solution_plotter(self,\
                         show_plots = False,\
                         collect_x_ghost_points = False,\
                         collect_y_ghost_points = False,\
                         **kwargs):
        """Append self.__python_tmp with a creation of the
        solution_plotter instance."""
        # Import the bout_plotters class
        self.__python_tmp +=\
            'from bout_runners.bout_plotters import solution_plotter\n'

        # Creates an instance of solution_and_error_plotter
        # Since we set qsub = True, the constructor will call the
        # collect_and_plot function
        # show_plots is set to false when running qsub
        self.__python_tmp +=\
          "make_my_plotter=solution_plotter(\n"+\
          "run_groups="             + str(self.__run_groups)        +",\n    "+\
          "show_plots="             + str(False)                  +",\n    "+\
          "collect_x_ghost_points=" + str(collect_x_ghost_points) +",\n    "+\
          "collect_y_ghost_points=" + str(collect_y_ghost_points) +",\n    "+\
          "directory='"             + str(self.__directory)        +"',\n    "+\
          "file_extension='"        + str(self.__file_extension)   +"',\n    "+\
          "variables="              + str(kwargs['variables'])    +",\n    "+\
          "plot_direction="         +\
                                str(kwargs['plot_direction'])     +",\n    "+\
          "plot_times="             + str(kwargs['plot_times'])   +",\n    "+\
          "number_of_overplots="    +\
          "qsub = True"                                          + ")\n"
#}}}
#}}}

#solution_and_error_plotter specific
#{{{
#{{{solution_and_error_plotter
    def solution_and_error_plotter(self,\
                         show_plots = False,\
                         collect_x_ghost_points = False,\
                         collect_y_ghost_points = False,\
                         **kwargs):
        """Append self.__python_tmp with a creation of the
        solution_and_error_plotter instance."""

        # Import the bout_plotters class
        self.__python_tmp +=\
            'from bout_runners.bout_plotters import solution_and_error_plotter\n'

        # Creates an instance of solution_and_error_plotter
        # Since we set qsub = True, the constructor will call the
        # collect_and_plot function
        # show_plots is set to false when running qsub
        self.__python_tmp +=\
          "make_my_plotter=solution_and_error_plotter(\n"+\
          "run_groups="             + str(self.__run_groups)        +",\n    "+\
          "show_plots="             + str(False)                  +",\n    "+\
          "collect_x_ghost_points=" + str(collect_x_ghost_points) +",\n    "+\
          "collect_y_ghost_points=" + str(collect_y_ghost_points) +",\n    "+\
          "directory='"             + str(self.__directory)        +"',\n    "+\
          "file_extension='"        + str(self.__file_extension)   +"',\n    "+\
          "variables="              + str(kwargs['variables'])    +",\n    "+\
          "plot_direction="         +\
                                str(kwargs['plot_direction'])     +",\n    "+\
          "plot_times="             + str(kwargs['plot_times'])   +",\n    "+\
          "number_of_overplots="    +\
                            str(kwargs['number_of_overplots'])    +",\n    "+\
          "qsub = True"                                          + ")\n"
#}}}
#}}}

#convergence_plot specific
#{{{
#{{{convergence_plotter
    def convergence_plotter(self,\
                         show_plots = False,\
                         collect_x_ghost_points = False,\
                         collect_y_ghost_points = False,\
                         **kwargs):
        """Append self.__python_tmp with a creation of the convergence_plotter
        instance."""

        # Import the bout_plotters class
        self.__python_tmp +=\
            'from bout_runners.bout_plotters import convergence_plotter\n'

        # Creates an instance of solution_and_error_plotter
        # Since we set qsub = True, the constructor will call the
        # collect_and_plot function
        # show_plots is set to false when running qsub
        self.__python_tmp +=\
          "make_my_plotter=convergence_plotter(\n"+\
          "run_groups="             + str(self.__run_groups)        +",\n    "+\
          "show_plots="             + str(False)                  +",\n    "+\
          "collect_x_ghost_points=" + str(collect_x_ghost_points) +",\n    "+\
          "collect_y_ghost_points=" + str(collect_y_ghost_points) +",\n    "+\
          "directory='"             + str(self.__directory)        +"',\n    "+\
          "file_extension='"        + str(self.__file_extension)   +"',\n    "+\
          "variables="              + str(kwargs['variables'])    +",\n    "+\
          "convergence_type='"      +\
                                str(kwargs['convergence_type'])  +"',\n    "+\
          "qsub = True"                                          + ")\n"
#}}}
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
