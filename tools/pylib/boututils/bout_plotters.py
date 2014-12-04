#!/usr/bin/env python
"""Classes and functions for plotting BOUT++ runs."""

# NOTE: This document uses folding. A hash-symbol followed by three {'s
# denotes the start of a fold, and a hash-symbol followed by three }'s
# denotes the end of a fold
__authors__ = 'Michael Loeiten'
__email__   = 'mmag@fysik.dtu.dk'
__version__ = '0.62beta'
__date__    = '27.11.2014'

import os
import re
from boututils import shell
from boututils.plot_style import set_style
from boututils.common_bout_functions import create_folder,\
                                            warning_printer,\
                                            check_for_plotters_errors,\
                                            wait_for_runs_to_finish,\
                                            find_variable_in_BOUT_inp
from boutdata import collect
import matplotlib.pyplot as plt
from pylab import plot
import numpy as np
import time
import re
import warnings

# TODO: Check if the l2 norm is correct if one takes into account that
#       the boundary point is staggered
# Consider: Make it easier to call manually
#           Make a demo for bout_plotters (idea: get_file_name to find
#           folder structure when called manually, use the same example as in)
# Consider: Remake get_time_string such that two legends cannot be the
#           same

#{{{class bout_plotter
# As an inherit class uses the super function, the class must allow an
# object as input
class bout_plotter(object):
    """Parent class for all plotter classes."""

#{{{__init__
    def __init__(self,\
        run_groups             = False,\
        show_plots             = False,\
        collect_x_ghost_points = False,\
        collect_y_ghost_points = False,\
        directory              = False,\
        file_extension         = 'png',\
        variables              = False,\
        qsub                   = False):
        """The constructor of this parent class is called by all plotter
        classes.
        
        If bout_plotters are called from a qsub routine, the variable 'qsub'
        will be set to 'True', which means that the collection and
        plotting routine will be called from this constructor"""

        print('\nNow making the plots\n')

        # Create an errorfilter to catch RuntimeWarings 
        # This will catch the division by zero Runtimewarning in np.log
        warnings.simplefilter("error", RuntimeWarning)
        # Member function which stores warnings
        self.warnings         = []

        # Common member functions for all plotters
        self.run_groups       = False
        # self.variables is a list of all the variables to be collected from
        # the dump file.
        self.variables              = variables
        self.show_plots             = show_plots
        self.collect_x_ghost_points = collect_x_ghost_points
        self.collect_y_ghost_points = collect_y_ghost_points
        self.directory              = directory
        self.file_extension         = file_extension

        # Assigning run_groups
        if qsub == False:
            # If the plotter are called manually it can be called in the
            # same manner as run_group is called from bout_runners.py, or as
            # a list of the folder one wishes to plot from

            # Therefore, we must convert self.run_groups so that no matter
            # how the run_groups are given, the run_groups will be on the
            # same form, that is
            # self.run_groups = {0:[dmp_folder_for_the_group], 1:...}

            # Rewrite manually made run_groups
            if type(run_groups) == list:
                self.run_groups = {}
                for group_no, job in enumerate(run_groups):
                    # Make the dictionary
                    self.run_groups[group_no] = [job]

            # If an instance have been made from bout_runners.py,
            # self.run_groups will have the following form
            # self.run_groups =\
            #   {0:{'dmp_folder':[...], 'job_status',[...]}, 1:...}

            # Rewrite run_groups made from bout_runners.py
            if type(run_groups[0]) == dict:                
                groups = run_groups.keys()
                self.run_groups = {}
                for group in groups:
                    # Make the dictionary
                    self.run_groups[group] =\
                        run_groups[group]['dmp_folder']
        else:
            # If the program has been called from the torque system,
            # qsub = True is given

            # If this is the case, collect_and_plot will be called from
            # the constructor (this is done to easen the script written
            # in the qsub runner in bout_runners.py)

            # Because of this, the reset_folder_counter_and_self_errors
            # in convergence_plotter is not called if we do convergence
            # runs from a cluster.
            if hasattr(self, 'convergence_type'):
                # Initialize the folder counter and self.errors
                self.reset_folder_counter_and_self_errors()

            # Wait for the runs to finish, and call collect_and_plot as
            # they finish
            wait_for_runs_to_finish(\
                run_groups, self.directory, self.group_done_function)
#}}}                

#{{{__del__ 
    def __del__(self):
        """The destructor will print all the error messages (if any)"""
        if len(self.warnings) == 0:
            print('\n'*3 + ' ' + '~'*70)
            print("| No WARNINGS detected before instance destruction in"+\
                  " 'bout_plotters'. |")
            print(' ' + '~'*70 + '\n'*3)
        else:
            print('\n'*3 + 'The following warnings were detected:')
            print('-'*80)
            for warning in self.warnings:
                print(warning + '\n')
            print('\n'*3)
#}}}        

# Functions called by the constructor
#{{{group_done_function
    def group_done_function(self, run_groups, group):
        """Makes a plot of a finished run group obtained from the input """
        # Rewrite self.run_groups to be the finished
        # group, and call collect_and_plot
        self.run_groups = {0:run_groups[group]['dmp_folder']}
        
        # Call the plotting routine
        self.collect_and_plot()
#}}}

# Common plotter functions
#{{{try_to_collect
    def try_to_collect(self, variable, path):
        """Tries to collect a variable. If an error occurs, it will be
        recasted to a warning, and the return value will be
        'skip_iteration'"""
        try:
            collected =\
                collect(variable, path=path, info=False,\
                        xguards = self.collect_x_ghost_points,\
                        yguards = self.collect_y_ghost_points)
        except ValueError:
            message = 'Could not collect from ' + path + '.'+\
                      ' Check that you are collecting the correct '+\
                      'variables. (Eventually also the .log (and/or '+\
                      '.err) file(s).'
            self.warnings.append(message)
            warning_printer(message)
            collected = 'skip_iteration'
        return collected
#}}}            
#}}}



#{{{class solution_plotter
class solution_plotter(bout_plotter):
    #{{{docstring
    """Takes a dictionary of the run_groups with the keys 'dmp_folder' and
    'status' as an input, collects the data and makes a solution plot for
    each job. Written so that it can be easily called from bout_runners.py

    The plots will be made with one column for the numerical solution.
    The number of rows equals the number of specified variables.
    """
    #}}}

# The constructor
#{{{__init__
    def __init__(self,\
        run_groups             = False,\
        show_plots             = False,\
        collect_x_ghost_points = False,\
        collect_y_ghost_points = False,\
        directory              = False,\
        file_extension         = 'png',\
        variables              = False,\
        plot_direction         = False,\
        plot_times             = False,\
        number_of_overplots    = False,\
        qsub                   = False):
        """Constructor for the solution_and_error_plotter"""

        # Declare child class' memberdata first, as the creator of
        # bout_plotter can run collect_and_plot()
        self.plot_direction      = plot_direction
        self.plot_times          = plot_times
        self.number_of_overplots = number_of_overplots

        # Data members which makes solution_plotter and
        # solution_and_error_plotter different
        try:
            # If this constructor is called from solution_and_error_plotter,
            # then the following member data is already set. In other words,
            # if self.plot_id is already set, we will skip setting the
            # following member data
            self.plot_id
        except AttributeError:
            # additional_subplot_index is effectively setting the number of
            # columns
            # Setting this to 20 gives two columns
            self.additional_subplot_index = 10
            self.style_name               = 'single_plot'
            self.plot_id                  = 'solution' 
        # Call the constructor of the superclass
        super(solution_plotter, self).__init__(\
                run_groups = run_groups,\
                show_plots = show_plots,\
                collect_x_ghost_points = collect_x_ghost_points,\
                collect_y_ghost_points = collect_y_ghost_points,\
                directory = directory,\
                file_extension = file_extension,\
                variables = variables,\
                qsub = qsub)

        # Check for errors
        plotter_error_checker =\
           check_for_plotters_errors('solution_plot',\
                                     variables = self.variables,\
                                     plot_direction = self.plot_direction,\
                                     plot_times = self.plot_times,\
                                     number_of_overplots = self.number_of_overplots) 
#}}}                

# Main function
#{{{collect_and_plot
    def collect_and_plot(self):
        """Drives the collection and plotting of the solution and error plot"""

        # Plotting preparation 
        number_of_variables = len(self.variables)

        # In this kind of plot, one group equals one job
        # We will here make one plot per job
        groups = self.run_groups.keys()
        for group in groups:
            # Each job is stored in a list
            # As we would like the job as a string we will pick the
            # zeroth element of the list
            job = self.run_groups[group][0]

            subplot_index = number_of_variables*100 +\
                            self.additional_subplot_index

            try:
                # Collect the time array, nout and timestep
                t = collect("t_array",\
                            xguards = self.collect_x_ghost_points,\
                            yguards = self.collect_y_ghost_points,\
                            path=job, info=False)
            except ValueError:
                message = "Could not collect 't_arry' in " + job + "."+\
                          " Check the log-file."
                self.warnings.append(message)
                warning_printer(message)
                # Jump to next iteration
                continue

            nout = t.shape[0]

            try:
                timestep = t[1]
            except IndexError:
                message = job + ' did not properly start.'
                self.warnings.append(message)
                warning_printer(message)
                # Jump to next iteration
                continue

            # Find how to slice
            x_slice, y_slice, z_slice = self.find_slices(job)
            # Jump to the next iteration if the collection in
            # find_slices failed
            if x_slice == 'skip_iteration':
                continue


            # Find t_array_indices
            t_array_indices =\
                self.get_t_array_indices(nout, timestep)

            # Call the style
            style = set_style(self.style_name, rows=len(self.variables))
            fig = plt.gcf()

            len_variables = len(self.variables)

            # Do the plotting
            for variable_no, variable in enumerate(self.variables):
                # Find the labels
                xlabel, ylabel = self.get_labels(x_slice, y_slice, z_slice) 

                if ylabel == '':
                    ylabel = variable

                solution = collect(variable,\
                                   xguards = self.collect_x_ghost_points,\
                                   yguards = self.collect_y_ghost_points,\
                                   path=job, info = False)
                # Collect additional variables (if any)
                additional_variables = self.collect_more(variable, job)


                subplot_index += 1
                for index in t_array_indices:
                    current_time = self.get_time_string(index, timestep)
                    ax = fig.add_subplot(subplot_index)
                    ax.plot(solution[index, x_slice, y_slice, z_slice],\
                             label="t="+current_time)
                if variable_no == 0:
                    ax.set_title('Numerical solution')
                ax.set_ylabel(ylabel)
                if variable_no + 1 == len_variables:
                    ax.set_xlabel(xlabel)
                style.plot_look_nice(ax)

                # If there are some additional variables
                if additional_variables != False:
                    # Plot the additional variables, and return the new
                    # subplot_index
                    subplot_index = self.additional_plot(\
                        t_array_indices, subplot_index, timestep,\
                        x_slice, y_slice, z_slice, len_variables,\
                        xlabel, additional_variables, variable_no, style, fig)

            variables = '_'.join(self.variables)
            plot_name = job + '/' + variables +\
                        '_' + self.plot_id + '_plot.' + self.file_extension

            plt.savefig(plot_name)
            print('Plot saved to ' + plot_name +'\n')

            # Move plot to comparison folder
            compare_folder = self.directory + '/compare_' + self.plot_id
            create_folder(compare_folder)
            dmp_folder_no_path = job.replace('/','-')
            command = 'cp ' + plot_name + ' ' + compare_folder +\
                      '/' + dmp_folder_no_path + '_' + variables +\
                      '.' + self.file_extension
            shell(command)
            if self.show_plots:
                plt.show()
#}}}

# Functions called by the main function
#{{{
#{{{collect_more
    def collect_more(self, *args):
        """Virtual function to be overridden by child classes."""
        return False
#}}}

#{{{additional_plot
    def additional_plot(self, **kwargs):
        """Virtual function to be overridden by child classes."""
        return False
#}}}        

#{{{find_slices
    def find_slices(self, in_folder):
        """Find how to slice the arrays obtained from 'collect' in order to make a plot.
        The arrays are indexed [t, x, y, z]

        self.plot_directions determines the slicing of the plot. Each
        direction can be set to a number in the grid, or to 'all'. Note that
        at least one direction has to be set to 'all', but no more than two.

        To be updated when 2D plots will be implimented.
        """
        
        # Used to find the lengths
        whole_array = self.try_to_collect(self.variables[0], in_folder)
        if whole_array == 'skip_iteration':
            # Jump to the next iteration
            return 'skip_iteration', None, None 

        x_len = len( whole_array[0,:,0,0] )
        y_len = len( whole_array[0,0,:,0] )
        z_len = len( whole_array[0,0,0,:] )

        number_of_slices = 0
        keys = self.plot_direction.keys()
        for key in keys:
            if key == 'x':
                if self.plot_direction[key] == 'all':
                    x_slice = slice(None)
                    number_of_slices += 1
                # Check if the index is valid
                elif self.plot_direction[key] <= x_len:
                    x_slice = self.plot_direction[key]
                else:
                    x_slice = x_len - 1
            if key == 'y':
                if self.plot_direction[key] == 'all':
                    y_slice = slice(None)
                    number_of_slices += 1
                # Check if the index is valid
                elif self.plot_direction[key] <= y_len:
                    y_slice = self.plot_direction[key]
                else:
                    y_slice = y_len - 1
            if key == 'z':
                if self.plot_direction[key] == 'all':
                    z_slice = slice(None)
                    number_of_slices += 1
                # Check if the index is valid
                elif self.plot_direction[key] <= x_len:
                    z_slice = self.plot_direction[key]
                else:
                    z_slice = z_len - 1
        
        return x_slice, y_slice, z_slice
#}}}

#{{{get_t_array_indices
    def get_t_array_indices(self, nout, timestep):
        """Get the indices of the over plots"""

        if self.number_of_overplots != False:
            # Divide into list of integers 
            step = int(np.floor(float(nout)/float(self.number_of_overplots)))
            # 1 is the smallest step in range
            if step == 0:
                step = 1
            t_array_indices = range( 0, nout, step )
            # Set the corresponding time values
            self.plot_times = [index*timestep for index in t_array_indices]
            # If get one more index than asked for
            while len(self.plot_times) > self.number_of_overplots:
                self.plot_times = self.plot_times[:-1]

        # Convert the time array to the corresponding indices
        t_array_indices = \
            self.convert_plot_times_to_indices(nout, timestep)

        return t_array_indices
#}}}

#{{{get_time_string
    def get_time_string(self, index, timestep):
        """Returns the time as an formatted string"""
        # Python string format cookbook:
        # http://mkaz.com/2012/10/10/python-string-format/
        current_time = index*timestep
        if (current_time > 100) or\
           ((current_time < 0.01) and (current_time != 0.0)):
            # Use scientific format
            return "{:.2e}".format(float(current_time))
        else:
            # Use standard format
            return "{:.2f}".format(float(current_time))
#}}}            

#{{{get_labels
    def get_labels(self, x_slice, y_slice, z_slice):
        """Obtain the x and y label for the plot"""
        ylabel = ''
        xlabel = ''
        if x_slice == slice(None):
            # The x grid should never be plotted on the y axis
            xlabel = 'x grid number'
        if y_slice == slice(None):
            if xlabel == '':
                xlabel = 'y grid number'
            else:
                ylabel = 'y grid number'
        if z_slice == slice(None):
            if xlabel == '':
                xlabel = 'z grid number'
            else:
                ylabel = 'z grid number'
        return xlabel, ylabel
#}}}
#}}}

# Auxiliary functions
#{{{
#{{{convert_plot_times_to_indices
    def convert_plot_times_to_indices(self, nout, timestep):
        """Convert self.plot_times into corresponding indices"""
        # Convert the plot_time to index in the t_array
        time_indices = []
        for the_time in self.plot_times:
            time_indices.append(round(the_time/timestep))
        for plot_times_index, t_array_index in enumerate(time_indices):
            # Check if there are numbers in the plot_times out of bounds
            if t_array_index >= nout:
                # If they are, they will be set to the last possible
                # index
                time_indices[plot_times_index] = nout - 1
        return time_indices
#}}}            
#}}}
#}}}



#{{{class solution_and_error_plotter
class solution_and_error_plotter(solution_plotter):
    #{{{docstring
    """Inherits from 'solution_plotter'

    The plots will be made with two columns, one for the numerical solution,
    and one for the errors.
    The number of rows equals the number of specified variables.

    NOTE: Requires that the runs have been run with mms=true.

    For more information, see 'solution_plotter'
    """
    #}}}

# The constructor
#{{{__init__
    def __init__(self,\
        run_groups             = False,\
        show_plots             = False,\
        collect_x_ghost_points = False,\
        collect_y_ghost_points = False,\
        directory              = False,\
        file_extension         = 'png',\
        variables              = False,\
        plot_direction         = False,\
        plot_times             = False,\
        number_of_overplots    = False,\
        qsub                   = False):
        """Constructor for the solution_and_error_plotter"""

        # Declare child class' memberdata first, as the creator of
        # bout_plotter can run collect_and_plot()
        # Setting this to 20 gives two columns
        self.additional_subplot_index = 20
        self.style_name               = 'two_columns'
        self.plot_id                  = 'solution_error' 

        # Call the constructor of the superclass
        super(solution_and_error_plotter, self).__init__(\
                run_groups = run_groups,\
                show_plots = show_plots,\
                collect_x_ghost_points = collect_x_ghost_points,\
                collect_y_ghost_points = collect_y_ghost_points,\
                directory = directory,\
                file_extension = file_extension,\
                variables = variables,\
                plot_direction = plot_direction,\
                plot_times = plot_times,\
                number_of_overplots = number_of_overplots,\
                qsub = qsub)
#}}}                

# Functions called by the main function
#{{{
#{{{collect_more
    def collect_more(self, variable, job):
        """Collect the error of the current variable."""
        error = collect("E_" + variable,\
                        xguards = self.collect_x_ghost_points,\
                        yguards = self.collect_y_ghost_points,\
                        path=job, info=False)
        return [error]
#}}}

#{{{additional_plot
    def additional_plot(self,\
                        t_array_indices, subplot_index, timestep,\
                        x_slice, y_slice, z_slice, len_variables,\
                        xlabel, additional_variables, variable_no, style, fig):
        """Plots the error."""
        subplot_index += 1
        error = additional_variables[0]
        for index in t_array_indices:
            current_time = self.get_time_string(index, timestep)
            ax = fig.add_subplot(subplot_index)
            ax.plot(error[index, x_slice, y_slice, z_slice],\
                     label="t="+current_time)
        if variable_no == 0:                         
            ax.set_title('Error')
        if variable_no + 1 == len_variables:
            ax.set_xlabel(xlabel)
            xlabel2_set = True
        style.plot_look_nice(ax)

        return subplot_index
#}}}
#}}}
#}}}



#{{{class convergence_plotter
class convergence_plotter(bout_plotter):
    #{{{docstring
    """Takes a dictionary of the run_groups with the keys 'dmp_folder' and
    'status' as an input, collects the data and makes a convergence plot for
    each of the run groups.
    Written so that it can be easily called from convergence_runner in
    bout_runners.py

    NOTE: Requires that the runs have been run with mms=true.
    """
    #}}}

# The constructor
#{{{__init__
    def __init__(self,\
        run_groups             = False,\
        show_plots             = False,\
        collect_x_ghost_points = False,\
        collect_y_ghost_points = False,\
        directory              = False,\
        file_extension         = 'png',\
        variables              = False,\
        convergence_type       = False,\
        qsub                   = False):
        """Constructor for the convergence_plotter """

        # Declare child class' memberdata first, as the creator of
        # bout_plotter can run collect_and_plot()
        self.convergence_type  = convergence_type
        self.style_name        = 'single_plot'

        # Call the constructor of the superclass
        # NOTE: If qsub = True, the parent class is going to call the
        # collecting and plotting routines
        super(convergence_plotter, self).__init__(\
            run_groups = run_groups,\
            show_plots = show_plots,\
            collect_x_ghost_points = collect_x_ghost_points,\
            collect_y_ghost_points = collect_y_ghost_points,\
            directory = directory,\
            file_extension = file_extension,\
            variables = variables,\
            qsub = qsub)

        # Initialize the folder counter and self.errors
        self.reset_folder_counter_and_self_errors()

        # If qsub=True, the check has already been made. Therefore, we
        # can place the super function above this as it is going to call
        # the collecting and plotting routines whenever qsub=True
        # Check for errors
        # Find the grid and timestep from the folder
        groups = run_groups.keys()
        for group in groups:
            # Find the folder name from run_groups
            folder_names = run_groups[group]['dmp_folder']
            # Appendable lists
            grids = []
            timestep = []
            for folder_name in folder_names:
                # Split the strings by /
                folder_list = folder_name.split('/')
                # Find the folder containing the mesh folder
                mesh_folder = [folder for folder in folder_list if\
                               'mesh' in folder or 'MZ' in folder]
                # As mesh_folder is stored as a list, we unpack the list
                # by taking the 0th element
                mesh_folder = mesh_folder[0]
                # Find the number in the mesh_folder
                # As findall returns a list, we unpack the list by
                # taking the 0th element
                # We are at the same time converting the number to a
                # integer number
                number = int(re.findall(r'\d+', mesh_folder)[0])
                # Add the number to the grid list
                grids.append(number)

                # Find the folder containing the timestep folder
                timestep_folder = [folder for folder in folder_list if\
                                   'timestep' in folder]
                # The folder may contain nout as well, so we remove this
                # If nout is in the folder, it is separated by a '_'
                # We take the zeroth element to get the string out of
                # the list
                timestep_folder = timestep_folder[0].split('_')
                # Select the timestep part of the folder
                timestep_folder = [folder for folder in timestep_folder if\
                                   'timestep' in folder]
                # The timestep number is separated by a '-'
                # We take the zeroth element to get the string out of
                # the list and remove 'timestep-' in order to select the number
                timestep_now = timestep_folder[0].replace('timestep-','')
                # Save the number as a float
                timestep.append(float(timestep_now))

            # Cast grids in to a dictionary (due to the construction of
            # the plotter_error_checker)
            grids = {'reconstructed_grids':grids}

            # Check for errors
            plotter_error_checker =\
               check_for_plotters_errors('convergence_plot',\
                                         variables = self.variables,\
                                         convergence_type =
                                         self.convergence_type,\
                                         timestep = timestep,\
                                         grids = grids)
#}}}                

# Main function
#{{{collect_and_plot
    def collect_and_plot(self):
        """Drives the collection and plotting of the convergence groups"""

        # Finding the convergence group from the dictionary
        groups = self.run_groups.keys()

        # Loop over the convergence groups
        for group_no in groups:
            folders_in_group = self.run_groups[group_no]
            no_of_folders_in_group = len(folders_in_group)

            # Loop over the folders in a convergence group
            for in_folder in folders_in_group:
                self.convergence_collect(in_folder)
                # If the grid_size and error values has been collected
                # for all folders in a group
                if self.folder_no_counter == no_of_folders_in_group:
                    # Make a plot
                    self.plot_and_print_convergence(in_folder)

                    # Reset the counter
                    self.reset_folder_counter_and_self_errors()
#}}}

# Functions called by the main function
#{{{
#{{{convergence_collect
    def convergence_collect(self, in_folder): 
        """Collects the data members belonging to a convergence plot"""

        self.folder_no_counter += 1
        # Collect the variables
        for variable in self.variables:
            # Collect if possible
            error_array = self.try_to_collect('E_'+variable, in_folder)
            if error_array == 'skip_iteration':
                # Jump to the next iteration
                continue

            # Check that there is a match between the desired and
            # obtained nout
            # Find nout
            # Check if nout is in the folder name
            if 'nout' in in_folder:
                # Split the string
                nout_desired = in_folder.split('_')
                # Get the nout
                nout_desired = [folder for folder in nout_desired\
                                if 'nout' in folder]
                # Find the numbers
                # Unpack the list by taking the zeroth element
                nout_desired = re.findall(r'\d+', nout_desired[0])
                # Convert the number to an integer
                nout_desired = int(nout_desired[0])
            else:
                # nout is found in BOUT.inp
                nout_desired = find_variable_in_BOUT_inp(in_folder, 'nout')
            # Find the obtained nout
            t = collect('t_array', path=in_folder,\
                        xguards = self.collect_x_ghost_points,\
                        yguards = self.collect_y_ghost_points,\
                        info=False)
            nout_obtained = len(t)
            # nout_obtained also includes the initial conditions
            if nout_desired != nout_obtained - 1:
                message = in_folder + ' did not finish all the timesteps.'
                self.warnings.append(message)
                warning_printer(message)
                # Jump to the next iteration
                continue
    
            # We have already found the E_inf norm for one time step. The
            # infinity error of this is the max of the absolute of the
            # error
            self.errors[variable]['error_2'].append(\
                np.sqrt(np.mean( error_array**2.0 )) \
                )
            self.errors[variable]['error_inf'].append(\
                np.max(np.abs( error_array )) \
                )
    
            # We want to find the spacing.
            folders = in_folder.split('/')
           
            if self.convergence_type == 'spatial':
                # The convergence run is constructed in a way such that all
                # of the gridspaces shoul be equal for a equidistant grid

                # Find the folder cotaining one of the grid spacing among
                # all the folders

                # Note that the guard cells are not collected as we are
                # only interested in the spacing
                for folder in folders:
                    if ('nx' in folder) or\
                       ('ny' in folder) or\
                       ('MZ' in folder):
                        if 'ny' in folder:
                            spacing =\
                                collect("dy", path=in_folder,
                                        xguards = False,\
                                        yguards = False,\
                                        info=False)
                            # We collect the maximum spacing
                            spacing = np.amax(spacing)
                        elif 'MZ' in folder:
                            spacing =\
                                collect("dz", path=in_folder,\
                                        xguards = False,\
                                        yguards = False,\
                                        info=False)
                                # There is only one possible grid spacing in
                                # MZ
                        else:
                            # nx is in folder
                            spacing =\
                                collect("dx", path=in_folder,\
                                        xguards = False,\
                                        yguards = False,\
                                        info=False)
                            # We collect the maximum spacing
                            spacing = np.amax(spacing)

                        self.errors[variable]['spacing'].append(spacing)
                        # Exit the for-loop
                        break
                        
            elif self.convergence_type == 'temporal':
                # Find the folder cotaining the timestep among all the folders
                for folder in folders:
                    # Since we are doing a temporal convergence test,
                    # timestep is bound to be in one of the folder names
                    if 'timestep' in folder:
                        # We want to extract the number
                        # The guard_time info is separated by a _
                        guard_time = folder.split('_')
                        for info in guard_time:
                            if 'timestep' in info:
                                # info is now separated like this
                                # timestep-<number>
                                timestep_number = info.split('-')
                                # The number is stored the least
                                number = float(timestep_number[-1])
                                self.errors[variable]['spacing'].append(number) 
                                # Exit the for-loop
                                break
    #}}}

#{{{plot_and_print_convergence
    def plot_and_print_convergence(self, in_folder):
        """Plots and prints the result from a convergence. Reset the
        variables belonging to the convergence when done"""

        # Loop over the given variables
        for variable in self.variables:

            # Obtaining the file_name
            file_name = self.get_convergence_path(variable, in_folder)
            # Sorting the data and obtain the order
            try:
                order_inf, order_2 = self.prepare_convergence_plot(variable, file_name)
                # Check that no RuntimeWarnings occured
                if order_inf == 'RuntimeWarning' or order_2 == 'RuntimeWarning':
                    continue
            except ValueError:
                message = "The order could not be found for '" +\
                          variable +"' belonging to " + in_folder + "." +\
                          " Check that at least two runs finished."
                self.warnings.append(message)
                warning_printer(message)
                # Jump to the next iteration
                continue
                        
            # Plot errors
            # Set the plotting style
            style = set_style(self.style_name)
            fig = plt.gcf()
            ax = fig.add_subplot(111)
            ax.plot(\
                     self.errors[variable]['spacing'],\
                     self.errors[variable]['error_inf'],\
                     'r-^',\
                     label=r'$L_\infty$')
            # In the log-log plot, we have 
            # ln(y) = a*ln(x) + ln(b)
            # y = error
            # x = spacing (found from linear regression)
            # a = order
            # b = where the straight line intersects with the ordinate
            #
            # Using the logarithmic identities
            # ln(x^a) = a*ln(x)
            # ln(xy) = ln(x) + ln(y)
            # ln(x/y) = ln(x) - ln(y)
            #
            # We usually find b for x = 0. Here, on the other hand, we find it
            # by solving the equation for the smallest grid point:
            # ln[y(x[-1])] = a*ln(x[-1]) + ln(b),
            # so
            # ln(b) = ln[y(x[-1])] - a*ln(x[-1])
            # =>
            # ln(y) = a*ln(x) - a*ln(x[-1]) + ln[y(x[-1])]
            #       = a*[ln(x)-ln(x[-1])] + ln[y(x[-1])]
            #       = a*[ln(x/x[-1])] + ln[ys(x[-1])]
            #       = ln[(x/x[-1])^a*y(x[-1])]
            ax.plot(\
                     (self.errors[variable]['spacing'][-1],\
                      self.errors[variable]['spacing'][0]),\
                     (\
                      ((self.errors[variable]['spacing'][-1]/\
                        self.errors[variable]['spacing'][-1])**order_inf)*\
                        self.errors[variable]['error_inf'][-1],\
                      ((self.errors[variable]['spacing'][0]/\
                        self.errors[variable]['spacing'][-1])**order_inf)*\
                        self.errors[variable]['error_inf'][-1]\
                     ),\
                     'm--',\
                     label=r"$\mathcal{O}_{L_\infty}"+"%.2f"%(order_inf)+r"$")
            ax.plot(\
                     self.errors[variable]['spacing'],\
                     self.errors[variable]['error_2'],\
                     'b-o',\
                     label=r'$L_2$')
            ax.plot(\
                     (self.errors[variable]['spacing'][-1],\
                      self.errors[variable]['spacing'][0]),\
                     (\
                      ((self.errors[variable]['spacing'][-1]/\
                        self.errors[variable]['spacing'][-1])**order_2)*\
                        self.errors[variable]['error_2'][-1],\
                      ((self.errors[variable]['spacing'][0]/\
                        self.errors[variable]['spacing'][-1])**order_2)*\
                        self.errors[variable]['error_2'][-1]\
                     ),\
                     'c--',\
                     label=r"$\mathcal{O}_{L_2}"+"%.2f"%(order_2)+r"$")
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlabel("Mesh spacing")
            ax.set_ylabel("Error norm")

            # Clean the plot
            style.plot_look_nice(ax)

            # Save the plot
            plt.savefig(file_name + '.' + self.file_extension)
            print('\nPlot saved to ' + file_name + '.' + self.file_extension +\
                  '\n'*2)
            # Move plot to comparison folder
            compare_folder = self.directory + '/compare_convergence'
            create_folder(compare_folder)
            file_name_no_path = file_name.replace('/','-')
            command = 'cp ' + file_name + '.' + self.file_extension +\
                      ' ' + compare_folder + '/' + file_name_no_path +\
                      '.' + self.file_extension
            shell(command)
            if self.show_plots:
                plt.show()
            plt.close()
#}}}
#}}}

# Auxiliary functions
#{{{
#{{{reset_folder_counter_and_self_errors
    def reset_folder_counter_and_self_errors(self):
        """Reset self.folder_no_counter and initializes the self.errors
        dictionary"""
        self.folder_no_counter = 0

        # Make a list from the list in self.variables
        dict_keys = [variable for variable in self.variables]
        # Initialize the dict
        self.errors = dict.fromkeys(dict_keys)
        # New set of dictionary keys
        dict_keys = ['error_2', 'error_inf', 'spacing']
        for variable in self.variables:
            # Initialize a dict for each of the variables
            self.errors[variable] = { key:[] for key in dict_keys}
#}}}

#{{{get_convergence_path
# Take the convergence type as input
    def get_convergence_path(self, variable, in_folder):
        """Finds the folder to put the convergence files in. Returns a
        relative file name"""
        # Make a list of the directories
        convergence_folders = in_folder.split('/')
        # Make a empty list of the folders we would like to remove
        # from the name
        rm = []
        # Empty list to fill what we remove
        guard_and_time = []
        grid           = []
        for folder in convergence_folders:
            if ('MXG' in folder) or\
               ('MYG' in folder) or\
               ('nout' in folder) or\
               ('timestep' in folder):
                rm.append(folder)
                guard_and_time.append(folder)
            elif ('MZ' in folder) or\
                 ('mesh-nx' in folder) or\
                 ('mesh-ny' in folder):
                rm.append(folder)
                grid.append(folder)

        # Remove the found folders from the convergence folder
        for element in rm:
            convergence_folders.remove(element)

        if self.convergence_type == 'spatial':
            # Add spatial_convergence to the folder structure
            convergence_folders.append('spatial_convergence')
            # The guard_and_time folder will be appended to the name
            if len(guard_and_time) != 0:
                convergence_name = guard_and_time[0] + '_'
            else:
                convergence_name = ''
        elif self.convergence_type == 'temporal':
            # Add temporal_convergence to the folder structure
            convergence_folders.append('temporal_convergence')
            # The grid folder will be appended to the name
            if len(grid) != 0:
                convergence_name = grid[0] + '_'
            else:
                convergence_name = ''

        # Convert the folder to a string
        convergence_folders = '/'.join(convergence_folders)
        # Create the folder
        create_folder(convergence_folders)

        # The final file name
        convergence_name = convergence_folders + '/' + \
                           convergence_name + variable
        return convergence_name
#}}}

#{{{prepare_convergence_plot
# Don't send in in_folder, but rather file_name
    def prepare_convergence_plot(self, variable, file_name):
        """Prepares the convergence plot by finding the filename and
        order."""
        # Sort the data in case it is unsorted (may be the case from
        # qsub runners)
        list_of_tuples_to_be_sorted =\
            zip(self.errors[variable]['spacing'],\
                self.errors[variable]['error_inf'],\
                self.errors[variable]['error_2'])
        # Sort the list
        # Note that we are sorting in reverse order, as we want the
        # highest grid spacing first
        sorted_list = sorted(list_of_tuples_to_be_sorted, reverse = True)
        # Unzip the sorted list
        self.errors[variable]['spacing'],\
            self.errors[variable]['error_inf'],\
            self.errors[variable]['error_2'] =\
            zip(*sorted_list)

        # Initialize the orders
        order_2 = [' '*7]
        order_inf = [' '*7]

        # The order will be found by finding a linear fit between two
        # nearby points in the error-spacing plot. Hence, we must let 
        # the index in the for loop run to the length minus one
        # If an runtimewarning is occuring
        for index in range(len(self.errors[variable]['spacing']) - 1):
            # p = polyfit(x,y,n) finds the coefficients of a polynomial p(x)
            # of degree that fits the data, p(x(i)) to y(i), in a least squares 
            # sense.
            # The result p is a row vector of length n+1 containing the
            # polynomial coefficients in descending powers

            # Check if the logarithm has bad values
            # Create a filter
            try:
                spacing_start   = np.log(self.errors[variable]['spacing'][index]) 
                spacing_end     = np.log(self.errors[variable]['spacing'][index + 1])
                error_start_2   = np.log(self.errors[variable]['error_2'][index])
                error_end_2     = np.log(self.errors[variable]['error_2'][index + 1])
                error_start_inf = np.log(self.errors[variable]['error_inf'][index])
                error_end_inf   = np.log(self.errors[variable]['error_inf'][index + 1])
            except RuntimeWarning:
                message = 'np.log returned a RuntimeWarning in ' + file_name + ', '+\
                          'when finding the logarithm of the error. ' +\
                          'Check if the error is zero.'
                self.warnings.append(message)
                warning_printer(message)
                return 'RuntimeWarning', 'RuntimeWarning'

            # Finding the order in the two norm
            order = np.polyfit([spacing_start, spacing_end],\
                               [error_start_2, error_end_2], 1)
            order_2.append(order[0])
            # Finding the infinity order
            order = np.polyfit([spacing_start, spacing_end],\
                               [error_start_inf, error_end_inf], 1)
            order_inf.append(order[0])

        outstring = zip(self.errors[variable]['spacing'],\
                        self.errors[variable]['error_2'],\
                        order_2,\
                        self.errors[variable]['error_inf'],\
                        order_inf)
        # Write the found orders
        f = open( file_name + '.txt', 'w' )        
        header = '#spacing    error_2    order_2    error_inf    order_inf' 
        # Write the header to a file and on screen
        f.write(header + '\n')
        print('\nNow printing the results of the convergence test:')
        print(header)
        for nr, line in enumerate(outstring):
            # Write the data to a file
            string = '    '.join(str(element) for element in line)
            f.write(string + "\n")
            # The first line contains a string which cannot be converted
            # to a float under the 'order_' sections
            if nr == 0:
                print("{:4.2e}    {:.2e}              {:.2e}".format(\
                  line[0], line[1], line[3]))
            else:
                print("{:.2e}    {:.2e}   {:.3e}  {:.2e}     {:.2e}".format(\
                      line[0], line[1], line[2], line[3], line[4]))
        f.close()
        print('\n')

        return order_inf[-1], order_2[-1]
#}}}
#}}}
#}}}
