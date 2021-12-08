def _resolve_slice(ind, num):
    if isinstance(ind, slice):
        ret = list(ind.indices(num))
        d = ret[1] - ret[0]
        rem = d % ret[2]
        if rem != 0:
            ret[1] -= rem
            ret[1] += ret[2]
        return ret
    else:
        if ind < 0:
            ind = num + ind
        if ind < 0 or ind >= num:
            raise IndexError("%d is out of range [%d:%d]" % (ind, 0, num))
        return [ind, ind + 1, 1]


def _resolve_slices(inds, nums):
    if len(inds) != len(nums):
        raise IndexError(
            "Expexted %d dimensions, but got %d dimensions." % (len(nums), len(inds))
        )
    ret = []
    for i in range(len(inds)):
        ret.append(_resolve_slice(inds[i], nums[i]))
    return ret
