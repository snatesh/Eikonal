import ctypes
import numpy as np
from multipledispatch import dispatch
from sys import exit

libngjquad = ctypes.CDLL('libngjquad.so')


def ngjQuad ( n, m, a, b, c, tol, tolc, 
              alph, use_newton, use_wolfe,
              nthreads ):
  gjquad = libngjquad.ngjquad_T2( n, m, a, b, c, tol, tolc,
                                  alph, use_newton, use_wolfe,
                                  nthreads )
  N = libngjquad.getN(gjquad) 
  Z0 = np.zeros((3*N,1), dtype=np.double, order='F')
  libngjquad.copyGJQuad ( gjquad, 
                          Z0.ctypes.data_as(
                            ctypes.POINTER(
                              ctypes.c_double)) )
  return Z0 

libngjquad.ngjquad_T2.argtypes = [ctypes.c_uint,\
                                  ctypes.c_uint,\
                                  ctypes.c_double,\
                                  ctypes.c_double,\
                                  ctypes.c_double,\
                                  ctypes.c_double,\
                                  ctypes.c_double,\
                                  ctypes.c_double,\
                                  ctypes.c_bool,\
                                  ctypes.c_bool,\
                                  ctypes.c_uint]
libngjquad.ngjquad_T2.restype = ctypes.c_void_p

libngjquad.copyGJQuad.argtypes = [ctypes.c_void_p,\
                                  ctypes.POINTER(ctypes.c_double)]

libngjquad.copyGJQuad.restype = None

libngjquad.getN.argtypes = [ctypes.c_void_p]
libngjquad.getN.restype =  ctypes.c_uint
