#!/usr/bin/env python
"""Functions common for files associated with bout_runners.py Could also
have interest independently of bout_runner associated files."""
from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object

# NOTE: This document uses folding. A hash-symbol followed by three {'s
# denotes the start of a fold, and a hash-symbol followed by three }'s
# denotes the end of a fold
__authors__ = 'Michael Loeiten'
__email__   = 'mmag@fysik.dtu.dk'
__version__ = '1.0beta'
__date__    = '26.02.2015'

import re
import os
import time
from boututils.shell import shell


#{{{find_variable_in_BOUT_inp
def find_variable_in_BOUT_inp(directory, variable):
    """Finds an integer value of a variable in BOUT.inp"""
    # Length of the variable we want to find
    variable_length = len(variable)
    
    # Open the file
    f = open(directory + '/BOUT.inp','r')
    # Read the file
    inp_file = f.read()
    # Close the file
    f.close()

    # Split the files to lines
    inp_lines = inp_file.split('\n')
    # Remove leading and ending spaces in the strings
    inp_lines = [line.strip() for line in inp_lines]
    # Seek for the match
    variable_line = [line for line in inp_lines if line[0:variable_length]
    == variable]
    
    # If find more than one line
    if len(variable_line) == 1:
        variable_line = variable_line[0]
        # Find the variable number
        variable_value = re.findall(r'\d+', variable_line)
        # Return the integer value
        return int(variable_value[0])
    elif len(variable_line) > 1:
        raise RuntimeError("Several lines starting with " + variable + \
                           " was found. Check your BOUT.inp file!")
    else:
        return 'Variable not found'
#}}}

#{{{move_log_err_files                            
def move_log_err_files(job, directory):
    """Move .err and .log output from the torque system"""
    # Move the error and log file
    command = 'mv ' + job + '.log ' + job + '.err ' +\
        directory + '/qsub_output/'
    shell(command)
#}}}

#{{{create_folder
def create_folder(folder):
    """Creates a folder if it doesn't exists"""
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(folder + " created\n")
#}}}

#{{{warning_printer
def warning_printer(message):
    """Printing warnings"""
    print('\n'*3 + '*'*37 + 'WARNING' + '*'*36)
    # Makes sure that no more than 80 characters are printed out at
    # the same time
    message_chunks=[]
    for chunk in message_chunker(message):
        rigth_padding = ' '*(76 - len(chunk))
        print('* ' + chunk + rigth_padding + ' *')
    print('*'*80 + '\n'*3)
#}}}

#{{{message_chunker
def message_chunker(message, chunk=76):
    """Generator used to chop a message so it doesn't exceed some
    width"""
    for start in range(0, len(message), chunk):
        yield message[start:start + chunk]
#}}}

#{{{wait_for_runs_to_finish
def wait_for_runs_to_finish(run_groups, directory,\
                            group_done_function = lambda *args, **kwargs:None,\
                            *args, **kwargs):
#{{{docstring
    """Checks if the run_groups submitted to the cluster are finished.
    
    If a job has finished, the .err and .log file will be moved to
    the qsub_output folder.
    
    Takes a function 'group_done_function' as an, which describes what to
    be done when a group has finished. If no function is given,
    nothing will be done (default is a lambda function which returns
    None)

    Terminates when all run_groups are done"""
#}}}

    # Checks whether or not all the run_groups has completed
    while len(list(run_groups.keys())) != 0:
        groups = list(run_groups.keys())

        # Checks whether or not a group has finished to run
        for group in groups:

            status = run_groups[group]['job_status']

            # NOTE: The following could be done nicely with all,
            # however, all(iterable) has a strange behavior in ipython

            # Count the number of dones by making a list
            # comprehension, and finding the length of the list
            a_list_of_done = [element for element in status\
                              if element=='done']
            no_of_done = len(a_list_of_done)

            # Chek if all runs in a group has finished
            if no_of_done == len(status):

                # If a function has been given, it will override the
                # default lambda function
                # Function which calls whatever needs to be done when a
                # group of jobs has finished
               
                group_done_function(run_groups = run_groups, group = group)

                # Remove the group from run_groups
                run_groups.pop(group)
            else:
                # If not all the runs in a convergence group has
                # finished, we will see if any new run_groups are finished
                # since the last time we checked
                # We check if a job has finished by checking if both
                # the .log file and .err file has been created from a
                # job

                # Get the files in the current folder as a list
                files = [f for f in os.listdir('.') if os.path.isfile(f)]
                for index_no, job in\
                    enumerate(\
                    run_groups[group]['job_status']):
                   
                    # If the job has finished, a .log and a .err is
                    # made
                    if (job + '.log' in files) and (job + '.err' in files):
                        # Move the qsub output
                        print(job + ' just finished\n')
                        move_log_err_files(job, directory)

                        # Set the job to done
                        run_groups[group][\
                            'job_status'][index_no] = 'done'

        # Sleep for five seconds
        time.sleep(5)
#}}}



#{{{class check_for_plotters_errors
class check_for_plotters_errors(object):
    """A class where the constructor searches for eventual user
    errors"""

#{{{Constructor
    def __init__(self, plot_type, **kwargs):
        """Constructor for the check_for_plotters_errors.
        The constructor calls the error checkers."""
        self.plot_type = plot_type

        if self.plot_type == 'solution_plot' or\
           self.plot_type == 'solution_and_error_plot':
            self.solution_plotter_error_checker(**kwargs)
        elif self.plot_type == 'convergence_plot':
            self.convergence_error_checker(**kwargs)
        else:
            raise TypeError ("The given 'plot_type' '" + str(self.plot_type) +\
                             "' is invalid. See run_with_plots"+\
                             " documentation for valid possibilities.")
#}}}

# Common functions for error checking of all plotters
#{{{check_keyword_arguments
    def check_keyword_arguments(self, keywords, options=False, **kwargs):
        """Check that the necessary additional keyword arguments are given"""

        kwarg_keys = list(kwargs.keys())
        kw_error = " needs to be given as a keyword argument when"+\
                   " running '" + self.plot_type + "'."

        # If options = 1, it means that this keyword must be given
        if options == 1:
            for keyword in keywords:
                if (keyword in kwarg_keys) == False:
                    raise RuntimeError ("'" + keyword + "'" + kw_error)

        # If options = 2, it means that either of the following keywords
        # must be given
        if options == 2:
            for keyword in keywords:
                error_message = "Either '" + keyword[0] + "' or '" +\
                                keyword[1] + "'" + kw_error
                if (keyword[0] in kwarg_keys)==False and\
                   (keyword[1] in kwarg_keys)==False:
                    raise RuntimeError (error_message)
        
#}}}

# Plotter specific error checkers
#{{{solution_plotter_error_checker
    def solution_plotter_error_checker(self, **kwargs):
        """Checks for specific solution_plotter and solution_and_error_plotter
        errors."""

        # Check that the necessary keyword arguments are given
        keywords = ['variables', 'plot_direction']
        self.check_keyword_arguments(keywords, options=1, **kwargs)


        if type(kwargs['variables']) != list:
            raise TypeError ("'variables' must be given as a list."+\
                             " See documentation")
        if type(kwargs['plot_direction']) != dict:
            raise TypeError ("'plot_direction' must be given as a dict."+\
                             " See documentation")

        keywords = [('plot_times', 'number_of_overplots')]
        self.check_keyword_arguments(keywords, options=2, **kwargs)
#}}}

#{{{convergence_error_checker
    def convergence_error_checker(self, **kwargs):
        """Checks for specific convergence_plotter errors"""

        # Check that the necessary keyword arguments are given
        keywords = ['variables', 'convergence_type']
        self.check_keyword_arguments(keywords, options=1, **kwargs)

        # Check that convergence is set correctly
        if (kwargs['convergence_type'] != 'temporal')\
            and (kwargs['convergence_type'] != 'spatial'):
            raise ValueError("The keyword 'convergence_type' must be set to"+\
                             " either 'spatial' or 'temporal'")

        if kwargs['convergence_type'] == 'spatial':            
            # Check if the conditions are set for a convergence test
            if kwargs['grids'] == False:
                raise TypeError("In order to make a convergence test, "+\
                                "you must specify a range of grid values.")
            keys = list(kwargs['grids'].keys())
            if len(kwargs['grids'][keys[0]]) <= 1:
                raise TypeError("In order to make a convergence test, "+\
                                "you must specify a range of grid values.")
        elif kwargs['convergence_type'] == 'temporal':
            # Check if a range is given, so that we can perform a
            # convergence test            
            if type(kwargs['timestep']) != list:
                raise TypeError("In order to make a convergence test, "+\
                                "you must specify a range of grid values.")
            if len(kwargs['timestep']) <= 1:
                raise TypeError("In order to make a convergence test, "+\
                                "you must specify a range of grid values.")
#}}}
#}}}
