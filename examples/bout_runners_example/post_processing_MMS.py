#!/usr/bin/env python

"""Post processing which performs MMS"""

from boutdata import collect
from boututils.showdata import showdata
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pylab import plot
import numpy as np
from _tkinter import TclError

#{{{perform_MMS_test
def perform_MMS_test(paths, extension='.pdf', show_plot=False):
    """Collects the data members belonging to a convergence plot"""

    # Make a variable to store the errors and the spacing
    data = {'error_2':[], 'error_inf':[], 'spacing':[]}

    # Loop over the runs in order to collect
    for path in paths:
        # Collect n_solution - n_numerical
        error_array = collect('E_n', path=path, info=False,\
                              xguards = False, yguards = False)
        # Pick the last time point
        error_array = error_array[-1]

        # The error in the 2-norm and infintiy-norm
        data['error_2']  .append( np.sqrt(np.mean( error_array**2.0 )) )
        data['error_inf'].append( np.max(np.abs( error_array )) )

        # Collect the spacings
        dx_spacing = collect("dx", path=path, info=False,\
                             xguards = False, yguards = False)
        dy_spacing = collect("dy", path=path, info=False,\
                             xguards = False, yguards = False)
        dz_spacing = collect("dz", path=path, info=False,\
                             xguards = False, yguards = False)
        # We are interested in the max of the spacing
        dx_spacing = np.max(dx_spacing)
        dy_spacing = np.max(dy_spacing)
        dz_spacing = np.max(dz_spacing)

        # TODO: Decide
        # Store the spacing in the data
        # data['spacing'].append(dx_spacing + dy_spacing + dz_spacing)
        # data['spacing'].append(dx_spacing)
        data['spacing'].append(dy_spacing)

    # Sort the data
    data = sort_data(data)

    # Find the order of convergence in the 2 norm and infinity norm
    order_2, order_inf = get_order(data)

    # Get the root name of the path (used for saving files)
    root_folder = paths[0].split('/')[0] + '/'

    # Get the name of the plot based on the first folder name
    name = 'MMS_'+path[0].replace('/', '_')

    # Print the convergence rate
    print_convergence_rate(data, order_2, order_inf, root_folder, name)

    # Plot
    # We want to show the lines of the last orders, so we send in
    # order_2[-1] and order_inf[-1]
    do_plot(data, order_2[-1], order_inf[-1],\
            root_folder, name, extension, show_plot)
#}}}

# Help functions
#{{{sort_data
def sort_data(data):
    """Sorts the data after highest grid spacing"""

    # Sort the data in case it is unsorted
    list_of_tuples_to_be_sorted =\
        list(zip(data['spacing'], data['error_inf'], data['error_2']))

    # Sort the list
    # Note that we are sorting in reverse order, as we want the
    # highest grid spacing first
    sorted_list = sorted(list_of_tuples_to_be_sorted, reverse = True)
    # Unzip the sorted list
    data['spacing'], data['error_inf'], data['error_2'] =\
            list(zip(*sorted_list))

    return data
#}}}

#{{{get_order
def get_order(data):
    # TODO: Check this
    # Initialize the orders
    order_2 = [' '*7]
    order_inf = [' '*7]

    # The order will be found by finding a linear fit between two
    # nearby points in the error-spacing plot. Hence, we must let
    # the index in the for loop run to the length minus one
    for index in range(len(data['spacing']) - 1):
        # p = polyfit(x,y,n) finds the coefficients of a polynomial p(x)
        # of degree that fits the data, p(x(i)) to y(i), in a least squares
        # sense.
        # The result p is a row vector of length n+1 containing the
        # polynomial coefficients in descending powers
        spacing_start   = np.log(data['spacing'][index])
        spacing_end     = np.log(data['spacing'][index + 1])
        error_start_2   = np.log(data['error_2'][index])
        error_end_2     = np.log(data['error_2'][index + 1])
        error_start_inf = np.log(data['error_inf'][index])
        error_end_inf   = np.log(data['error_inf'][index + 1])
        # Finding the order in the two norm
        order = np.polyfit([spacing_start, spacing_end],\
                           [error_start_2, error_end_2], 1)
        # Append it to the order_2
        order_2.append(order[0])

        # Finding the infinity order
        order = np.polyfit([spacing_start, spacing_end],\
                           [error_start_inf, error_end_inf], 1)
        # Append it to the order_inf
        order_inf.append(order[0])

    return order_2, order_inf
#}}}

#{{{print_convergence_rate
def print_convergence_rate(data, order_2, order_inf, root_folder, name):
    "Prints the convergence rates to the screen and to a file"
    outstring = list(zip(data['spacing'],\
                         data['error_2'],\
                         order_2,\
                         data['error_inf'],\
                         order_inf))

    # Write the found orders
    f = open( name + '.txt', 'w' )
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
#}}}

#{{{do_plot
def do_plot(data, order_2, order_inf, root_folder, name, extension, show_plot):
    """Function which handles the actual plotting"""

    # Plot errors
    # Set the plotting style
    title_size = 30
    plt.rc("font", size = 30)
    plt.rc("axes", labelsize = 25, titlesize = title_size)
    plt.rc("xtick", labelsize = 25)
    plt.rc("ytick", labelsize = 25)
    plt.rc("legend", fontsize = 30)
    plt.rc("lines", linewidth = 3.0)
    plt.rc("lines", markersize = 20.0)
    plt_size = (10, 7)
    fig_no = 1
    # Try to make a figure with the current backend
    try:
        fig = plt.figure(fig_no, figsize = plt_size)
    except TclError:
        # Switch if a backend needs the display
        plt.switch_backend('Agg')
        fig = plt.figure(fig_no, figsize = plt_size)

    ax = fig.add_subplot(111)

    # Plot errore
    # Plot the error-space plot for the 2-norm
    ax.plot(data['spacing'], data['error_2'], 'b-o', label=r'$L_2$')
    # Plot the error-space plot for the inf-norm
    ax.plot(data['spacing'], data['error_inf'], 'r-^', label=r'$L_\infty$')

    # Plot the order
    #{{{ Explanaition of the calculation
    # In the log-log plot, we have
    # ln(y) = a*ln(x) + ln(b)
    # y = error
    # x = spacing (found from linear regression)
    # a = order
    # b = where the straight line intersects with the ordinate
    #
    # Using the logarithmic identities
    # ln(x^a) = a*ln(x)
    # ln(x*y) = ln(x) + ln(y)
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
    #}}}
    # Order in the inf norm
    ax.plot(\
             (data['spacing'][-1],\
              data['spacing'][0]),\
             (\
              ((data['spacing'][-1] / data['spacing'][-1])**order_inf)*\
                data['error_inf'][-1],\
              ((data['spacing'][0]  / data['spacing'][-1])**order_inf)*\
                data['error_inf'][-1]\
             ),\
             'm--',\
             label=r"$\mathcal{O}_{L_\infty}"+"%.2f"%(order_inf)+r"$")
    # Order in the 2 norm
    ax.plot(\
             (data['spacing'][-1],\
              data['spacing'][0]),\
             (\
              ((data['spacing'][-1] / data['spacing'][-1])**order_2)*\
                data['error_2'][-1],\
              ((data['spacing'][0]  / data['spacing'][-1])**order_2)*\
                data['error_2'][-1]\
             ),\
             'c--',\
             label=r"$\mathcal{O}_{L_2}"+"%.2f"%(order_2)+r"$")

    # Set logaraithmic scale
    ax.set_yscale('log')
    ax.set_xscale('log')

    # Set axis label
    ax.set_xlabel("Mesh spacing")
    ax.set_ylabel("Error norm")

    # Make the plot look nice
    # Plot the legend
    leg = ax.legend(loc="best", fancybox = True, numpoints=1)
    leg.get_frame().set_alpha(0.5)
    # Plot the grid
    ax.grid()
    # Make sure no collision between the ticks
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    # Includes the xlabel if outside
    plt.tight_layout()

    # Save the plot
    plt.savefig(root_folder + name  + '.' + extension)
    print('\nPlot saved to ' + name + '.' + extension + '\n'*2)

    plt.show()
    plt.close()
#}}}
