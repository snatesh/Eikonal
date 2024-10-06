import ctypes
import numpy as np
from multipledispatch import dispatch
from sys import exit

libsfactors = ctypes.CDLL('libsfactors.so')

@dispatch(int, float, float)
def sFactors(N, a, b):
  """
  Structural constants for Jacobi polynomials on [-1,1]
  """
 
  H = np.zeros((N,1), dtype=np.double, order='F')
  libsfactors.sFactors_T1 ( N, a, b, 
                            H.ctypes.data_as(
                              ctypes.POINTER(
                                ctypes.c_double)) )
  return H 


@dispatch(int, float, float, float)
def sFactors(N, a, b, c):
  """
  Structural constants for Jacobi polynomials on
  the standard triangle
  """

  H = np.zeros((N,N), dtype=np.double, order='F')
  libsfactors.sFactors_T2 ( N, a, b, c, 
                            H.ctypes.data_as(
                              ctypes.POINTER(
                                ctypes.c_double)) )
  return H


@dispatch(int, int, int, float, float, float, float) 
def sFactors(n, k, j, a, b, c, d):
  """
  Structural constants for Jacobi polynomials on
  the standard tetrahedron
  """

  if n < k or k < j or j < 0: 
    exit("sFactors: Range Error (n >= k >= j >= 0)")
  return libsfactors.sFactors_T3(n, k, j, a, b, c, d) 


######################################################################

libsfactors.sFactors_T1.argtypes = [ctypes.c_uint,\
                                    ctypes.c_double,\
                                    ctypes.c_double,\
                                    ctypes.POINTER(ctypes.c_double)] 
libsfactors.sFactors_T1.restype = None


libsfactors.sFactors_T2.argtypes = [ctypes.c_uint,\
                                   ctypes.c_double,\
                                   ctypes.c_double,\
                                   ctypes.c_double,\
                                   ctypes.POINTER(ctypes.c_double)] 
libsfactors.sFactors_T2.restype = None

libsfactors.sFactors_T3.argtypes = [ctypes.c_uint,\
                                   ctypes.c_uint,\
                                   ctypes.c_uint,\
                                   ctypes.c_double,\
                                   ctypes.c_double,\
                                   ctypes.c_double,\
                                   ctypes.c_double] 
libsfactors.sFactors_T3.restype = ctypes.c_double


