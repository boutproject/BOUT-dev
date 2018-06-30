"""Wrapper for ndarray with extra attributes for BOUT++ fields.

"""

import numpy


class BoutArray(numpy.ndarray):
    """Wrapper for ndarray with extra attributes for BOUT++ fields.

    Parameters
    ----------
    input_array : array_like
        Data to convert to BoutArray
    attributes : dict
        Dictionary of extra attributes for BOUT++ fields

        Notably, these attributes should contain
        ``bout_type``. Possible values are:

        - scalar
        - Field2D
        - Field3D

        If the variable is an evolving variable (i.e. has a time
        dimension), then it is appended with a "_t"

    """

    # See https://docs.scipy.org/doc/numpy-1.13.0/user/basics.subclassing.html
    # for explanation of the structure of this numpy.ndarray wrapper

    def __new__(cls, input_array, attributes={}):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = numpy.asarray(input_array).view(cls)
        # add the dict of attributes to the created instance
        obj.attributes = attributes
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # ``self`` is a new object resulting from
        # ndarray.__new__(BoutArray, ...), therefore it only has
        # attributes that the ndarray.__new__ constructor gave it -
        # i.e. those of a standard ndarray.
        #
        # We could have got to the ndarray.__new__ call in 3 ways:
        # From an explicit constructor - e.g. BoutArray():
        #    obj is None
        #    (we're in the middle of the BoutArray.__new__
        #    constructor, and self.attributes will be set when we return to
        #    BoutArray.__new__)
        if obj is None:
            return
        # From view casting - e.g arr.view(BoutArray):
        #    obj is arr
        #    (type(obj) can be BoutArray)
        # From new-from-template - e.g boutarray[:3]
        #    type(obj) is BoutArray
        #
        # Note that it is here, rather than in the __new__ method, that we set
        # the default value for 'attributes', because this method sees all
        # creation of default objects - with the BoutArray.__new__ constructor,
        # but also with arr.view(BoutArray).
        self.attributes = getattr(obj, 'attributes', None)
        # We do not need to return anything

    def __format__(self, str):
        try:
            return super().__format__(str)
        except TypeError:
            return float(self).__format__(str)
