# cython: profile=True
# filename: host_props.pyx
# cython: boundscheck=False, overflowcheck=False

import numpy as np
import scipy.sparse

cimport numpy as np
cimport libc.stdlib as stdlib
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def mean_group_prop(gal_prop, gal_group_ID, group_ID):
    
    #get the arguments into shape
    cdef np.ndarray[np.float64_t, ndim=1] cgal_prop = np.ascontiguousarray(gal_prop,dtype=np.float64)
    cdef np.ndarray[np.intp_t, ndim=1] cgal_group_ID = np.ascontiguousarray(gal_group_ID,dtype=np.int)
    cdef np.ndarray[np.intp_t, ndim=1] cgroup_ID = np.ascontiguousarray(group_ID,dtype=np.int)
    
    #define some useful things
    cdef int N = len(gal_group_ID)
    cdef int i,j
    
    #the results
    Ngroup = len(group_ID)
    cdef np.ndarray[np.float64_t, ndim=1] result = np.zeros((Ngroup,), dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] m = np.zeros((Ngroup,), dtype=np.float64)

    for i in range(0,N):
        result[gal_group_ID[i]] += gal_prop[i]
        m[gal_group_ID[i]] += 1.0
    
    result = result/m    
    
    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def sum_group_prop(gal_prop, gal_group_ID, group_ID):
    pass



@cython.boundscheck(False)
@cython.wraparound(False)
def dist_group_prop(gal_prop, bins, gal_group_ID, group_ID):
    pass