import ctypes
import numpy as np
from multipledispatch import dispatch
from sys import exit
from sFactors import * 
from ngjQuad import *  
from pMat import *
from legendre import *
from jPoly import *
import matplotlib.pyplot as plt
libeikonal = ctypes.CDLL('libeikonal.so')

class eikonal(object):
  """ 
  Modal/spectral solver for Eikonal equation on Triangle
  We assume homogoeneous dirichlet boundary conditions
  """ 
  def __init__( self, _rhs, _u0, _a = 0.5, _b = 0.5, _c = 0.5,
                _wghted=False, _n = 4, _m = 6, _nthreads = 6,
                _qfname=None, _r = 0.005, _nthetas = 20 ):
    
    if _n <= 0:  
      exit("eikonal : Range Error ( n > 1) ")
    if _a <= -0.5 or _b <= -0.5 or _c <= -0.5:
      exit("eikonal : Range Error (a,b,c > -1/2)")

    ## jacobi parameters / method order / thread config
    self.a = _a
    self.b = _b
    self.c = _c
    self.weighted = _wghted
    self.n = _n
    self.m = _m
    self.nthreads = _nthreads
    self.N = int(0.5 * (_n) * (_n + 1))
    self.qfname = _qfname

    ## load or generate quadrature rule
    if not hasfile(self.qfname):
      self.qfname = 'triquadleg_n' + str(_n-1) + '_m' \
                     + str(_m-1) + '_N' + str(self.N) + '.txt'
      self.Z0 = ngjQuad ( _n, _m, _a, _b, _c, 1e-14, 1e-14, 
                          0, 0, 0, _nthreads )
      np.savetxt(self.qfname, self.Z0)
    else:
      self.Z0 = np.loadtxt(self.qfname)
    
    self.qN = int(len(self.Z0) / 3.0)
    self.X = self.Z0[0:self.qN]
    self.Y = self.Z0[self.N:2*self.qN]
    self.W = self.Z0[2*self.N:3*self.qN]
    # solution variables
    self.cu = self.uCoeffs(_u0)
    self.rhs = _rhs   
    # boundary nodes
    self.nl = self.N
    self.nb = self.N
    self.nh = self.N
    self.xl = np.zeros((self.nl,))
    self.yl, _ = legendre(self.nl); self.yl = (self.yl + 1.0) / 2.0 
    self.xb, _ = legendre(self.nb); self.xb = (self.xb + 1.0) / 2.0 
    self.yb = np.zeros((self.nb,))
    self.xh, _ = legendre(self.nh); self.xh = (self.xh + 1.0) / 2.0 
    self.yh = 1.0 - self.xh
    self.Xe = np.concatenate((self.xl, self.xb, self.xh))
    self.Ye = np.concatenate((self.yl, self.yb, self.yh))  
    self.Ne = self.nl + self.nb + self.nh

 
    ## derivative operators
    if not self.weighted: 
      self.Habc = sFactors(_n+1, _a, _b, _c)  
      self.Ha1bc = sFactors(_n+1, _a+1, _b, _c)
      self.Ha1b1c = sFactors(_n+1, _a+1, _b+1, _c)
      self.Ha1b1c1 = sFactors(_n+1, _a+1, _b+1, _c+1)
      self.Ha1bc1 = sFactors(_n+1, _a+1, _b, _c+1)
      self.Hab1c = sFactors(_n+1, _a, _b+1, _c);
      self.Hab1c1 = sFactors(_n+1, _a, _b+1, _c+1)
    else:
      self.Habc = sFactors(_n+2, _a, _b, _c)  
      self.Ha1bc = sFactors(_n+2, _a+1, _b, _c)
      self.Ha1b1c = sFactors(_n+2, _a+1, _b+1, _c)
      self.Ha1b1c1 = sFactors(_n+2, _a+1, _b+1, _c+1)
      self.Ha1bc1 = sFactors(_n+2, _a+1, _b, _c+1)
      self.Hab1c = sFactors(_n+2, _a, _b+1, _c);
      self.Hab1c1 = sFactors(_n+2, _a, _b+1, _c+1)

    self.Kabc_a1bc = kMat(_a, _b, _c, self.Habc, self.Ha1bc, _n-1, 0)
    self.Ka1bc_a1b1c = kMat(_a+1, _b, _c, self.Ha1bc, self.Ha1b1c, _n-1, 1)
    self.Ka1b1c_a1b1c1 = kMat(_a+1, _b+1, _c, self.Ha1b1c, self.Ha1b1c1, _n-1, 2)
    # promotion for RHS
    self.K = np.asfortranarray(
              self.Ka1b1c_a1b1c1.dot(self.Ka1bc_a1b1c.dot(self.Kabc_a1bc)))

    # promotion so derivs are in same basis
    self.K_a1bc1_a1b1c1 = np.asfortranarray(
                            kMat(_a+1, _b, _c+1, self.Ha1bc1, self.Ha1b1c1, _n-1, 1))
    self.K_ab1c1_a1b1c1 = np.asfortranarray(
                            kMat(_a, _b+1, _c+1, self.Hab1c1, self.Ha1b1c1, _n-1, 0))
    
    if not self.weighted:
      self.Dx = np.asfortranarray(
                self.K_a1bc1_a1b1c1.dot(
                  dMat  ( self.a, self.b, self.c,
                          self.Habc, self.Ha1bc1, _n-1, 0 ) ) )
      self.Dy = np.asfortranarray(
                self.K_ab1c1_a1b1c1.dot(
                  dMat  ( self.a, self.b, self.c,
                          self.Habc, self.Hab1c1, _n-1, 1 ) ) )
      self.polyE = jPoly(self.Ne, _n, _a, _b, _c, _nthreads, self.weighted)
      self.poly = jPoly(self.N, _n, _a, _b, _c, _nthreads, self.weighted)
      self.polypt = jPoly(1, _n, _a, _b, _c, 1, self.weighted)
      self.polyCirc = jPoly(_nthetas, _n, _a, _b, _c, _nthreads, self.weighted)
      self.polyl = jPoly(self.nl, _n, _a, _b, _c, 1, self.weighted)
      self.polyb = jPoly(self.nb, _n, _a, _b, _c, 1, self.weighted)
      self.polyh = jPoly(self.nh, _n, _a, _b, _c, 1, self.weighted)
    else: 
      self.Lx = lMat(self.a, self.b+1, self.c, self.Hab1c, self.Habc, _n, 1)
      self.Ly = lMat(self.a+1, self.b, self.c, self.Ha1bc, self.Habc, _n, 0)
      self.Wx = dMat(self.a+1, self.b+1, self.c+1, self.Ha1b1c1, self.Hab1c, _n, 0, True)
      self.Wy = dMat(self.a+1, self.b+1, self.c+1, self.Ha1b1c1, self.Ha1bc, _n, 1, True)
      self.Dx = np.asfortranarray(self.Lx.dot(self.Wx))
      self.Dy = np.asfortranarray(self.Ly.dot(self.Wy))
      self.polyE = jPoly(self.Ne, _n, _a+1, _b+1, _c+1, _nthreads, self.weighted)
      self.poly = jPoly(self.N, _n, _a+1, _b+1, _c+1, _nthreads, self.weighted)
      self.polypt = jPoly(1, _n, _a+1, _b+1, _c+1, 1, self.weighted)
      self.polyCirc = jPoly(_nthetas, _n, _a+1, _b+1, _c+1, _nthreads, self.weighted)
      self.polyl = jPoly(self.nl, _n, _a+1, _b+1, _c+1, 1, self.weighted)
      self.polyb = jPoly(self.nb, _n, _a+1, _b+1, _c+1, 1, self.weighted)
      self.polyh = jPoly(self.nh, _n, _a+1, _b+1, _c+1, 1, self.weighted)
    
    
    self.rhscoeffs = np.asfortranarray( self.rhsCoeffs() )
    self.G = np.transpose(self.Dx).dot(self.Dx) + np.transpose(self.Dy).dot(self.Dy)
    self.G = np.asfortranarray(self.G[0:self.N,0:self.N])
    
    ## boundary evaluator
    self.Ve = np.asfortranarray(computeV(self.polyE, self.N, self.Xe, self.Ye)) 

    ## post / path processing
    self.V = self.computeV(self.X, self.Y)
    self.nthetas = _nthetas
    self.thetas = np.linspace(0, 2*np.pi, _nthetas) 
    self.r = _r
    self.Xcirc = np.asfortranarray(_r * np.cos(self.thetas))
    self.Ycirc = np.asfortranarray(_r * np.sin(self.thetas))

  def solveP(self):
    libeikonal.eikonalSolveP(self.N,
                             self.Ne,
                             self.Dx.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)), 
                             self.Dy.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)),
                             self.G.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)),
                             self.rhscoeffs.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)), 
                             self.Ve.ctypes.data_as(
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


  def rhsCoeffs(self):
    cF, _ = computeCoeffs ( self.rhs, self.n-1, self.a, 
                            self.b, self.c, self.X, self.Y, 
                            self.W, self.nthreads)
    return cF


  def distToHyp(self,X,Y):
    Np = np.size(X,0)
    dh = np.zeros_like(X)
    for j in range(Np):
      x = np.array([X[j],Y[j]])
      z = np.array([(X[j]-Y[j] + 1) /2, 
                    (Y[j]-X[j] + 1) / 2])
      dh[j] = np.linalg.norm(x - z)
    return dh
    
  
  def distToTri(self,X,Y):
    Z = self.distToHyp(X,Y)
    mins = np.zeros_like(X)
    N = np.size(X,0)
    for j in range(N):
      mins[j] = np.min([X[j],Y[j],Z[j]])
    return mins
  
  def uCoeffs(self, u0):
    cu, _ = computeCoeffs ( u0, self.n-1, self.a, 
                            self.b, self.c, self.X, self.Y, 
                            self.W, self.nthreads )
    return cu


def hasfile(fn):
  try:
    open(fn, "r")
  except IOError:
    print("Generating " + fn)
    return 0
  print("Reading quad rule from " + fn)
  return 1

libeikonal.eikonalSolveH.argtypes = [ctypes.c_uint,\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double)]
libeikonal.eikonalSolveH.restype = None                                    

libeikonal.eikonalSolveP.argtypes = [ctypes.c_uint,\
                                     ctypes.c_uint,\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double)]
libeikonal.eikonalSolveP.restype = None     

libeikonal.getCoeffs.argtypes = [ctypes.c_void_p,\
                                 ctypes.POINTER(ctypes.c_double)]
libeikonal.getCoeffs.restype = None

libeikonal.createSolver.argtypes = [ctypes.c_uint,\
                                    ctypes.c_uint,\                               
                                    ctypes.c_uint,\                               
                                    ctypes.c_uint,\                               
                                    ctypes.c_uint,\                               
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.POINTER(ctypes.c_double),\
                                    ctypes.c_uint]
libeikonal.createSolver.restype = ctypes.c_void_p

