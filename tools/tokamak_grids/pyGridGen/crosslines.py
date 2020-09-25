from __future__ import print_function
from __future__ import division
from past.utils import old_div
import numpy as np
from numpy.lib.stride_tricks import as_strided
import itertools

def unique(a, atol=1e-08):
    """Find unique rows in 2d array

    Parameters
    ----------
    a : 2d ndarray, float
        array to find unique rows in
    atol : float, optional
        tolerance to check uniqueness

    Returns
    -------
    out : 2d ndarray, float
        array of unique values

    Notes
    -----
    Adapted to include tolerance from code at https://stackoverflow.com/questions/8560440/removing-duplicate-columns-and-rows-from-a-numpy-2d-array#answer-8564438

    """

    if np.issubdtype(a.dtype, float):
        order = np.lexsort(a.T)
        a = a[order]
        diff = np.diff(a, axis=0)
        np.abs(diff,out = diff)
        ui = np.ones(len(a), 'bool')
        ui[1:] = (diff > atol).any(axis=1)
        return a[ui]
    else:
        order = np.lexsort(a.T)
        a = a[order]
        diff = np.diff(a, axis=0)
        ui = np.ones(len(a), 'bool')
        ui[1:] = (diff != 0).any(axis=1)
        return a[ui]

def linelineintersect(a, b, atol=1e-08):
    """Find all intersection points of two lines defined by series of x,y pairs

    Intersection points are unordered
    Colinear lines that overlap intersect at any end points that fall within the overlap

    Parameters
    ----------
    a, b : ndarray
        2 column ndarray of x,y values defining a two dimensional line.  1st
        column is x values, 2nd column is x values.

    Notes
    -----
    An example of 3 segment line is: a = numpy.array([[0,0],[5.0,3.0],[4.0,1]])
    Function faster when there are no overlapping line segments


    add some lines for preventing zero-division
    """

    def x_y_from_line(line):
        """1st and 2nd column of array"""
        return (line[:, 0],line[:, 1])
    def meshgrid_as_strided(x, y, mask=None):
        """numpy.meshgrid without copying data (using as_strided)"""
        if mask is None:
            return (as_strided(x, strides=(0, x.strides[0]), shape=(y.size, x.size)),
                    as_strided(y, strides=(y.strides[0],0), shape=(y.size, x.size)))
        else:
            return (np.ma.array(as_strided(x, strides=(0, x.strides[0]), shape=(y.size, x.size)), mask=mask),
                    np.ma.array(as_strided(y, strides=(y.strides[0],0), shape=(y.size, x.size)), mask=mask))

    #In the following the indices i, j represents the pairing of the ith segment of b and the jth segment of a
    #e.g. if ignore[i,j]==True then the ith segment of b and the jth segment of a cannot intersect
    ignore = np.zeros([b.shape[0]-1, a.shape[0]-1], dtype=bool)

    x11, x21 = meshgrid_as_strided(a[:-1, 0], b[:-1, 0], mask=ignore)
    x12, x22 = meshgrid_as_strided(a[1: , 0], b[1: , 0], mask=ignore)
    y11, y21 = meshgrid_as_strided(a[:-1, 1], b[:-1, 1], mask=ignore)
    y12, y22 = meshgrid_as_strided(a[1: , 1], b[1: , 1], mask=ignore)

    #ignore segements with non-overlappiong bounding boxes
    ignore[np.ma.maximum(x11, x12) < np.ma.minimum(x21, x22)] = True
    ignore[np.ma.minimum(x11, x12) > np.ma.maximum(x21, x22)] = True
    ignore[np.ma.maximum(y11, y12) < np.ma.minimum(y21, y22)] = True
    ignore[np.ma.minimum(y11, y12) > np.ma.maximum(y21, y22)] = True

    #find intersection of segments, ignoring impossible line segment pairs when new info becomes available
    denom_ = np.empty(ignore.shape, dtype=float)
    denom = np.ma.array(denom_, mask=ignore)
    denom_[:, :] = ((y22 - y21) * (x12 - x11)) - ((x22 - x21) * (y12 - y11))
    #denom_tmp  = ((y22 - y21) * (x12 - x11)) - ((x22 - x21) * (y12 - y11)) # H.SETO
    denom_[:, :] = np.where(denom_ == 0.0, 1.e-100,denom_)

    ua_ = np.empty(ignore.shape, dtype=float)
    ua = np.ma.array(ua_, mask=ignore)
    ua_[:, :] = (((x22 - x21) * (y11 - y21)) - ((y22 - y21) * (x11 - x21)))
    ua_[:, :] /= denom

    ignore[ua < 0] = True
    ignore[ua > 1] = True

    ub_ = np.empty(ignore.shape, dtype=float)
    ub = np.ma.array(ub_, mask=ignore)
    ub_[:, :] = (((x12 - x11) * (y11 - y21)) - ((y12 - y11) * (x11 - x21)))
    ub_[:, :] /= denom


    ignore[ub < 0] = True
    ignore[ub > 1] = True

    nans_ = np.zeros(ignore.shape, dtype = bool)
    nans = np.ma.array(nans_, mask = ignore)
    nans_[:,:] = np.isnan(ua)

    if not np.ma.any(nans):

        #remove duplicate cases where intersection happens on an endpoint
#        ignore[np.ma.where((ua[:, :-1] == 1) & (ua[:, 1:] == 0))] = True
#        ignore[np.ma.where((ub[:-1, :] == 1) & (ub[1:, :] == 0))] = True
        ignore[np.ma.where((ua[:, :-1] < 1.0 + atol) & (ua[:, :-1] > 1.0 - atol) & (ua[:, 1:] < atol) & (ua[:, 1:] > -atol))] = True
        ignore[np.ma.where((ub[:-1, :] < 1 + atol) & (ub[:-1, :] > 1 - atol) & (ub[1:, :] < atol) & (ub[1:, :] > -atol))] = True

        xi = x11 + ua * (x12 - x11)
        yi = y11 + ua * (y12 - y11)
        return np.ma.compressed(xi), np.ma.compressed(yi)
    else:
        n_nans = np.ma.sum(nans)
        n_standard = np.ma.count(x11) - n_nans
        #I initially tried using a set to get unique points but had problems with floating point equality

        #create n by 2 array to hold all possible intersection points, check later for uniqueness
        points = np.empty([n_standard + 2 * n_nans, 2], dtype = float) #each colinear segment pair has two intersections, hence the extra n_colinear points

        #add standard intersection points
        xi = x11 + ua * (x12 - x11)
        yi = y11 + ua * (y12 - y11)
        points[:n_standard,0] = np.ma.compressed(xi[~nans])
        points[:n_standard,1] = np.ma.compressed(yi[~nans])
        ignore[~nans]=True


        #now add the appropriate intersection points for the colinear overlapping segments
        #colinear overlapping segments intersect on the two internal points of the four points describing a straight line
        #ua and ub have already serverd their purpose. Reuse them here to mean:
            #ua is relative position of 1st b segment point along segment a
            #ub is relative position of 2nd b segment point along segment a
        use_x = x12 != x11 #avoid vertical lines diviiding by zero
        use_y = ~use_x

        ua[use_x] = old_div((x21[use_x]-x11[use_x]), (x12[use_x]-x11[use_x]))
        ua[use_y] = old_div((y21[use_y]-y11[use_y]), (y12[use_y]-y11[use_y]))

        ub[use_x] = old_div((x22[use_x]-x11[use_x]), (x12[use_x]-x11[use_x]))
        ub[use_y] = old_div((y22[use_y]-y11[use_y]), (y12[use_y]-y11[use_y]))

        np.ma.clip(ua, a_min=0,a_max = 1, out = ua)
        np.ma.clip(ub, a_min=0,a_max = 1, out = ub)

        xi = x11 + ua * (x12 - x11)
        yi = y11 + ua * (y12 - y11)
        points[n_standard:n_standard + n_nans,0] = np.ma.compressed(xi)
        points[n_standard:n_standard + n_nans,1] = np.ma.compressed(yi)

        xi = x11 + ub * (x12 - x11)
        yi = y11 + ub * (y12 - y11)
        points[n_standard+n_nans:,0] = np.ma.compressed(xi)
        points[n_standard+n_nans:,1] = np.ma.compressed(yi)

        points_unique = unique(points)
        return points_unique[:,0], points_unique[:,1]


def find_inter( contour1, contour2):
    xi = np.array([])
    yi = np.array([])

    i=0
    ncombos = (sum([len(x.get_paths()) for x in contour1.collections]) *
            sum([len(x.get_paths()) for x in contour2.collections]))
    for linecol1, linecol2 in itertools.product(contour1.collections, contour2.collections):
        for path1, path2 in itertools.product(linecol1.get_paths(),linecol2.get_paths()):
            i += 1
            print('line combo %5d of %5d' % (i, ncombos))
            xinter, yinter = linelineintersect(path1.vertices, path2.vertices)

            xi = np.append(xi, xinter)
            yi = np.append(yi, yinter)

    #plt.scatter(xi, yi, s=20)
    #plt.show()
    return xi, yi
