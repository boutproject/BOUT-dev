def resolve_slice(ind,num):
    if isinstance(ind,slice):
        ret=list(ind.indices(num))
        d=ret[1]-ret[0]
        rem=d%ret[2]
        if rem != 0:
            ret[1]-=rem
            ret[1]+=ret[2]
        return ret
    else:
        if ind < 0:
            ind=num+ind
        if ind < 0 or ind >= num:
            raise IndexError()
        return [ind,ind+1,1]
        
def resolve_slices(inds,nums):
    if len(inds)!= len(nums):
        raise IndexError()
    ret=[]
    for i in range(len(inds)):
        ret.append(resolve_slice(inds[i],nums[i]))
    return ret
