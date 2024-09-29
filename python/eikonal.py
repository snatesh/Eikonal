import ctypes
import numpy as np
from multipledispatch import dispatch
from sys import exit
from sFactors import * 
from ngjQuad import *  

libdmat = ctypes.CDLL('libdmat.so')
libkmat = ctypes.CDLL('libkmat.so')

class eikonal(object):
  """ 
  Operators for derivatives in x,y 
  acting on coefficients of an order n Jacobi
  Polynomial expansion (on the standard triangle)
  """ 
  #def __init__(self, _a, _b, _c, _n = 10, _m = 16, _nthreads = 6):
  def __init__(self, _a, _b, _c, _n = 4, _m = 6, _nthreads = 6):
    
    if _n <= 0:  
      exit("eikonal : Range Error ( n > 1) ")
    if _a <= -0.5 or _b <= -0.5 or _c <= -0.5:
      exit("eikonal : Range Error (a,b,c > -1/2)")

    # jacobi parameters
    self.a = _a
    self.b = _b
    self.c = _c
    self.n = _n
    self.m = _m
    self.nthreads = _nthreads
    self.N = int(0.5 * (_n + 1) * (_n + 2))
    self.fname = 'triquadleg_n' + str(_n-1) + '_m' \
                  + str(_m-1) + '_N' + str(self.N) + '.txt'
    if not hasfile(self.fname):
      self.Z0 = ngjQuad ( _n, _m, _a, _b, _c, 1e-14, 1e-14, 
                          0, 0, 0, _nthreads )
      np.savetxt(self.fname, self.Z0)
    else:
      Z0 = np.loadtxt(self.fname)
    

    self.Habc = sFactors(_n+1, _a, _b, _c)  
    self.Ha1bc = sFactors(_n+1, _a+1, _b, _c)
    self.Ha1b1c = sFactors(_n+1, _a+1, _b+1, _c)
    self.Ha1b1c1 = sFactors(_n+1, _a+1, _b+1, _c+1)
    self.Ha1bc1 = sFactors(_n+1, _a+1, _b, _c+1)
    self.Hab1c1 = sFactors(_n+1, _a, _b+1, _c+1)

    self.Kabc_a1bc = kMat(_a, _b, _c, self.Habc, self.Ha1bc, _n-1, 0)
    self.Ka1bc_a1b1c = kMat(_a+1, _b, _c, self.Ha1bc, self.Ha1b1c, _n-1, 1)
    self.Ka1b1c_a1b1c1 = kMat(_a+1, _b+1, _c, self.Ha1b1c, self.Ha1b1c1, _n-1, 2)
    self.K = self.Ka1b1c_a1b1c1.dot(self.Ka1bc_a1b1c.dot(self.Kabc_a1bc))    

    self.Dx = self.K.dot(
                dMat  ( self.a, self.b, self.c,
                        self.Habc, self.Ha1bc1, _n-1, 0 ) )

    self.Dy = self.K.dot(
                dMat  ( self.a, self.b, self.c,
                        self.Habc, self.Hab1c1, _n-1, 1 ) )

    print(self.Dx)
    print(self.Dy)
    np.savetxt("Dx.txt",self.Dx)
    np.savetxt("Dy.txt",self.Dy)

def dMat(a, b, c, H, H1, n, mode):
  
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
                    ctypes.c_double)) )
  return D 

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
  return K 

def hasfile(fn):
  try:
    open(fn, "r")
  except IOError:
    print("Generating " + fn)
    return 0
  print("Reading quad rule from " + fn)
  return 1

libdmat.dMat.argtypes = [ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.c_uint,\
                         ctypes.c_uint,\
                         ctypes.POINTER(ctypes.c_double)] 
libdmat.dMat.restype = None

libkmat.kMat.argtypes = [ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.c_double,\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.POINTER(ctypes.c_double),\
                         ctypes.c_uint,\
                         ctypes.c_uint,\
                         ctypes.POINTER(ctypes.c_double)] 
libkmat.kMat.restype = None
