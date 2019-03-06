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
    if len(errors) != len(grid_spacing):
        raise ValueError("errors (len: {}) and grid_spacing (len: {}) should be the same length"
                         .format(len(errors), len(grid_spacing)))

    full_range = polyfit(log(grid_spacing), log(errors), 1)

    small_spacing = log(errors[-2] / errors[-1]) / log(grid_spacing[-2] / grid_spacing[-1])

    return (full_range[0], small_spacing)


def check_order(error_list, expected_order, tolerance=2.e-1, spacing=None):
    """Check if the actual_order is sufficiently close to the
    expected_order within a given tolerance

    """

    if len(error_list) < 2:
        raise RuntimeError("Expected at least 2 data points to calculate error")

    success=True
    for i in range(len(error_list)-1):
        if spacing is None:
            actual_order = log(errors[i] / errors[i+1]) / log(2)
        else:
            actual_order = log(errors[i] / errors[i+1]) / log(grid_spacing[i] / grid_spacing[i+1])

        if not isclose(actual_order, expected_order, atol=tolerance, rtol=0):
            success=False
    return success


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
    if len(errors) != len(grid_sizes):
        raise ValueError("errors (len: {}) and grid_sizes (len: {}) should be the same length"
                         .format(len(errors), len(grid_sizes)))

    dx = 1. / array(grid_sizes)
    message = "{}:\nGrid points | Error    | Rate\n".format(label)
    for i, grid_size in enumerate(grid_sizes):
        message += "{:<11} | {:f} | ".format(grid_size, errors[i])
        if i > 0:
            message += "{:f} \n".format(log(errors[i] / errors[i-1]) / log(dx[i] / dx[i-1]))
        else:
            message += "--\n"
    return message
