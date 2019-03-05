"""Functions for checking the error scaling of MMS or MES results

"""

from numpy import array, isclose, log, polyfit


def get_order(grid_spacing, errors):
    """Get the convergence order of errors over the full range of
    grid_spacing, and at small spacings

    Parameters
    ----------
    grid_spacing : list of float
        The grid spacing or inverse of number of grid points
    errors : list of float
        The error at each grid spacing

    Returns
    -------
    tuple of float
        The first value is the error scaling over the full range of
        grid spacings; the second value is the scaling over the last
        two points

    """
    full_range = polyfit(log(grid_spacing), log(errors), 1)

    small_spacing = log(errors[-2] / errors[-1]) / log(grid_spacing[-2] / grid_spacing[-1])

    return (full_range[0], small_spacing)


def check_order(actual_order, expected_order, tolerance=2.e-1):
    """Check if the actual_order is sufficiently close to the
    expected_order within a given tolerance

    """
    return isclose(actual_order, expected_order, rtol=tolerance)


def error_rate_table(errors, grid_sizes, label):
    """Create a nicely formatted table of the error convergence rate over
    the grid_sizes

    The error rate is calculated between adjacent points

    Parameters
    ----------
    errors : list of float
        The errors at each grid size
    grid_sizes : list of int
        The number of grid points
    label : string
        What the error is measuring

    Returns
    -------
    string

    """
    dx = 1. / array(grid_sizes)
    message = "{}:\nGrid points | Error    | Rate\n".format(label)
    for i in range(len(grid_sizes)):
        message += "{:<11} | {:f} | ".format(grid_sizes[i], errors[i])
        if i > 0:
            message += "{:f} \n".format(log(errors[i] / errors[i-1]) / log(dx[i] / dx[i-1]))
        else:
            message += "--\n"
    return message
