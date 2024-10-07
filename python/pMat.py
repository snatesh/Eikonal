import ctypes
import numpy as np


libdmat = ctypes.CDLL('libdmat.so')
libkmat = ctypes.CDLL('libkmat.so')
liblmat = ctypes.CDLL('liblmat.so')

def dMat(a, b, c, H, H1, n, mode, weighted=False):
  
  N = int(0.5 * (n + 1) * (n + 2))
  D = np.zeros((N,N), dtype=np.double, order='F')
  libdmat.dMat  ( a, b, c,  
                  H.ctypes.data_as(
                    ctypes.POINTER(
                    ctypes.c_double)),
                  H1.ctypes.data_as(
                    ctypes.POINTER(
                    ctypes.c_double)),
                  n, mode, 
                  D.ctypes.data_as(
                    ctypes.POINTER(
                    ctypes.c_double)),
                  weighted )
  return np.asfortranarray(D) 

def lMat(a, b, c, H, H1, n, mode):
  
  N = int(0.5 * (n + 1) * (n + 2))
  L = np.zeros((N,N), dtype=np.double, order='F')
  liblmat.lMat  ( a, b, c,  
                  H.ctypes.data_as(
                    ctypes.POINTER(
                    ctypes.c_double)),
                  H1.ctypes.data_as(
                    ctypes.POINTER(
                    ctypes.c_double)),
                  n, mode, 
                  L.ctypes.data_as(
                    ctypes.POINTER(
                    ctypes.c_double)) )
  return np.asfortranarray(L) 

def kMat(a, b, c, H, H1, n, mode):
  
  N = int(0.5 * (n + 1) * (n + 2))
  K = np.zeros((N,N), dtype=np.double, order='F')
  libkmat.kMat  ( a, b, c,  
                  H.ctypes.data_as(
                    ctypes.POINTER(
                    ctypes.c_double)),
                  H1.ctypes.data_as(
                    ctypes.POINTER(
                    ctypes.c_double)),
                  n, mode, 
                  K.ctypes.data_as(
                    ctypes.POINTER(
                    ctypes.c_double)) )
  return np.asfortranarray(K) 

libdmat.dMat.argtypes = [ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.c_uint,\
                         ctypes.c_uint,\
                         ctypes.POINTER(ctypes.c_double),
                         ctypes.c_bool] 
libdmat.dMat.restype = None

liblmat.lMat.argtypes = [ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.c_uint,\
                         ctypes.c_uint,\
                         ctypes.POINTER(ctypes.c_double)]

liblmat.lMat.restype = None

libkmat.kMat.argtypes = [ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.c_uint,\
                         ctypes.c_uint,\
                         ctypes.POINTER(ctypes.c_double)] 
libkmat.kMat.restype = None
