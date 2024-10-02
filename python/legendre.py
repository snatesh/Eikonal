import ctypes
import numpy as np



liblegquad = ctypes.CDLL('liblegquad.so')

def legendre(n):
  legx = np.zeros((n,))
  legw = np.zeros((n,))
  quad = liblegquad.legendre(n)
  liblegquad.copyQuad(quad, 
                      legx.ctypes.data_as(
                        ctypes.POINTER(
                        ctypes.c_double)),
                      legw.ctypes.data_as(
                        ctypes.POINTER(
                        ctypes.c_double)))
  liblegquad.deleteQuad(quad)
  return legx, legw 



liblegquad.legendre.argtypes = [ctypes.c_uint]
liblegquad.legendre.restype = ctypes.c_void_p

liblegquad.copyQuad.argtypes = [ctypes.c_void_p,\
                                ctypes.POINTER(ctypes.c_double),\
                                ctypes.POINTER(ctypes.c_double)]
liblegquad.copyQuad.restype = None

liblegquad.deleteQuad.argtypes = [ctypes.c_void_p]
liblegquad.deleteQuad.restype = None
