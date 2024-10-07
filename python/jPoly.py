import ctypes
import numpy as np

libjpoly = ctypes.CDLL('libjpoly.so')


def computeCoeffs(f, n, a, b, c, X, Y ,W, nthreads):
  Npts = np.size(X,0)
  N = int(0.5 * (n+1) * (n+2))
  poly = jPoly(Npts, n, a, b, c, nthreads)
  V = computeV(poly, N, X, Y)
  return np.transpose(V).dot(f(X,Y) * W), V
  

def computeV(poly, N, x, y):
  Npts = np.size(x,0)
  V = np.zeros((Npts, N), order='F')
  libjpoly.computeV(  poly, 
                      x.ctypes.data_as(
                        ctypes.POINTER(
                        ctypes.c_double)),
                      y.ctypes.data_as(
                        ctypes.POINTER(
                        ctypes.c_double)),
                      V.ctypes.data_as(
                        ctypes.POINTER(
                        ctypes.c_double)) )
  return np.asfortranarray(V)

def jPoly(Nx, n, a, b, c, nthreads, weighted=False):
  poly = libjpoly.jPoly_T2(Nx, n, a, b, c, nthreads, weighted)
  return poly


libjpoly.jPoly_T2.argtypes = [ctypes.c_uint,\
                              ctypes.c_uint,\
                              ctypes.c_double,\
                              ctypes.c_double,\
                              ctypes.c_double,\
                              ctypes.c_uint,
                              ctypes.c_bool]
libjpoly.jPoly_T2.restype = ctypes.c_void_p

libjpoly.computeV.argtypes = [ctypes.c_void_p,\
                              ctypes.POINTER(ctypes.c_double),\
                              ctypes.POINTER(ctypes.c_double),\
                              ctypes.POINTER(ctypes.c_double)]
libjpoly.computeV.restype = None
