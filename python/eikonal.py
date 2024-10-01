import nlopt
import ctypes
import numpy as np
from multipledispatch import dispatch
from sys import exit
from sFactors import * 
from ngjQuad import *  

libdmat = ctypes.CDLL('libdmat.so')
libkmat = ctypes.CDLL('libkmat.so')
libjpoly = ctypes.CDLL('libjpoly.so')
libeikonal = ctypes.CDLL('libeikonal.so')

class eikonal(object):
  """ 
  Operators for derivatives in x,y 
  acting on coefficients of an order n Jacobi
  Polynomial expansion (on the standard triangle)
  """ 
  def __init__(self, _a, _b, _c, _n = 4, _m = 6, _nthreads = 6, _r = 0.005, _nthetas = 100):
    
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
    self.N = int(0.5 * (_n) * (_n + 1))
    self.cu = np.asfortranarray(np.zeros((self.N,)))
    self.cu[0] = 1.0; # initialize to constant poly

    self.fname = 'triquadleg_n' + str(_n-1) + '_m' \
                  + str(_m-1) + '_N' + str(self.N) + '.txt'
    if not hasfile(self.fname):
      self.Z0 = ngjQuad ( _n, _m, _a, _b, _c, 1e-14, 1e-14, 
                          0, 0, 0, _nthreads )
      np.savetxt(self.fname, self.Z0)
    else:
      self.Z0 = np.loadtxt(self.fname)
    
    self.X = self.Z0[0:self.N]
    self.Y = self.Z0[self.N:2*self.N]
    self.W = self.Z0[2*self.N:3*self.N]
    
    self.poly = jPoly(self.N, _n, _a, _b, _c, _nthreads)
    self.V = self.computeV(self.X, self.Y)
    
    self.Habc = sFactors(_n+1, _a, _b, _c)  
    self.Ha1bc = sFactors(_n+1, _a+1, _b, _c)
    self.Ha1b1c = sFactors(_n+1, _a+1, _b+1, _c)
    self.Ha1b1c1 = sFactors(_n+1, _a+1, _b+1, _c+1)
    self.Ha1bc1 = sFactors(_n+1, _a+1, _b, _c+1)
    self.Hab1c1 = sFactors(_n+1, _a, _b+1, _c+1)

    self.Kabc_a1bc = kMat(_a, _b, _c, self.Habc, self.Ha1bc, _n-1, 0)
    self.Ka1bc_a1b1c = kMat(_a+1, _b, _c, self.Ha1bc, self.Ha1b1c, _n-1, 1)
    self.Ka1b1c_a1b1c1 = kMat(_a+1, _b+1, _c, self.Ha1b1c, self.Ha1b1c1, _n-1, 2)
    # promotion for RHS
    self.K = self.Ka1b1c_a1b1c1.dot(self.Ka1bc_a1b1c.dot(self.Kabc_a1bc))    
    # promotion so derivs are in same basis
    self.K_a1bc1_a1b1c1 = kMat(_a+1, _b, _c+1, self.Ha1bc1, self.Ha1b1c1, _n-1, 1)
    self.K_ab1c1_a1b1c1 = kMat(_a, _b+1, _c+1, self.Hab1c1, self.Ha1b1c1, _n-1, 0)

    self.Dx = np.asfortranarray(
              self.K_a1bc1_a1b1c1.dot(
                dMat  ( self.a, self.b, self.c,
                        self.Habc, self.Ha1bc1, _n-1, 0 ) ) )
    self.Dy = np.asfortranarray(
              self.K_ab1c1_a1b1c1.dot(
                dMat  ( self.a, self.b, self.c,
                        self.Habc, self.Hab1c1, _n-1, 1 ) ) )

    self.Cf = self.finvsq()
    self.edge_zeros = np.zeros_like(self.X)
    self.vl = self.computeV(self.edge_zeros, self.Y)
    self.vb = self.computeV(self.X, self.edge_zeros)
    self.vh = self.computeV(self.X, 1.0 - self.X)
    self.polypt = jPoly(1, _n, _a, _b, _c, 1)
    self.nthetas = _nthetas
    self.polyCirc = jPoly(_nthetas, _n, _a, _b, _c, _nthreads)
    self.thetas = np.linspace(0, 2*np.pi, _nthetas) 
    self.r = _r
    self.Xcirc = _r * np.cos(self.thetas)
    self.Ycirc = _r * np.sin(self.thetas) 

    print(np.isfortran(self.vl))

  def solve(self):
    libeikonal.eikonalSolve(self.N, 
                            self.Dx.ctypes.data_as(
                              ctypes.POINTER(
                              ctypes.c_double)), 
                            self.Dy.ctypes.data_as(
                              ctypes.POINTER(
                              ctypes.c_double)), 
                            self.vl.ctypes.data_as(
                              ctypes.POINTER(
                              ctypes.c_double)), 
                            self.vb.ctypes.data_as(
                              ctypes.POINTER(
                              ctypes.c_double)), 
                            self.vh.ctypes.data_as(
                              ctypes.POINTER(
                              ctypes.c_double)), 
                            self.cu.ctypes.data_as(
                              ctypes.POINTER(
                              ctypes.c_double)))
    return self.cu 

    
 
  def evalPoly(self, xpt, ypt):
    x = np.zeros((1,))
    y = np.zeros((1,))
    x[0] = xpt
    y[0] = ypt
    V = np.zeros((1, self.N), order='F')
    libjpoly.computeV(  self.polypt, 
                        x.ctypes.data_as(
                          ctypes.POINTER(
                          ctypes.c_double)),
                        y.ctypes.data_as(
                          ctypes.POINTER(
                          ctypes.c_double)),
                        V.ctypes.data_as(
                          ctypes.POINTER(
                          ctypes.c_double)) )
    return V

  def computeV(self, x, y):
    return computeV(self.poly, self.N, x, y)


  # constant speed in medium  
  # coefficient representation is just first  
  # standard basis vector
  def finvsq(self):
    cf = np.zeros((self.N,))
    cf[0] = 1.0
    return np.outer(cf,cf)
   
  def F(self, cu):
    vecx = self.Dx.dot(cu)
    vecy = self.Dy.dot(cu)
    return np.linalg.norm(np.outer(vecx,vecx) + 
                          np.outer(vecy,vecy) - 
                          self.Cf, ord='fro')**2
 

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
  return V

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

def jPoly(Nx, n, a, b, c, nthreads):
  poly = libjpoly.jPoly_T2(Nx, n, a, b, c, nthreads)
  return poly



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
  
libjpoly.jPoly_T2.argtypes = [ctypes.c_uint,\
                              ctypes.c_uint,\
                              ctypes.c_double,\
                              ctypes.c_double,\
                              ctypes.c_double,\
                              ctypes.c_uint]
libjpoly.jPoly_T2.restype = ctypes.c_void_p

libjpoly.computeV.argtypes = [ctypes.c_void_p,\
                              ctypes.POINTER(ctypes.c_double),\
                              ctypes.POINTER(ctypes.c_double),\
                              ctypes.POINTER(ctypes.c_double)]
libjpoly.computeV.restype = None

libeikonal.eikonalSolve.argtypes = [ctypes.c_uint,\
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.POINTER(ctypes.c_double)]
libeikonal.eikonalSolve.restype = None                                    
