import ctypes
import numpy as np
from multipledispatch import dispatch
from sys import exit
from sFactors import * 
from ngjQuad import *  
from pMat import *
from legendre import *
from jPoly import *

libeikonal = ctypes.CDLL('libeikonal.so')

class eikonal(object):
  """ 
  Operators for derivatives in x,y 
  acting on coefficients of an order n Jacobi
  Polynomial expansion (on the standard triangle)
  """ 
  def __init__( self, _rhs, _a = 0.5, _b = 0.5, _c = 0.5, 
                _n = 4, _m = 6, _nthreads = 6, 
                _r = 0.001, _nthetas = 50 ):
    
    if _n <= 0:  
      exit("eikonal : Range Error ( n > 1) ")
    if _a <= -0.5 or _b <= -0.5 or _c <= -0.5:
      exit("eikonal : Range Error (a,b,c > -1/2)")

    ## jacobi parameters / method order / thread config
    self.a = _a
    self.b = _b
    self.c = _c
    self.n = _n
    self.m = _m
    self.nthreads = _nthreads
    self.N = int(0.5 * (_n) * (_n + 1))
    # solution variables
    self.cu = self.cuCoeffs()#np.zeros((self.N,))
    # initialize to something parabaloid-y with y=0 plane
    #self.cu[0] = 1.0;
    #self.cu[1] = 0.0;
    #self.cu[2] = 1.0;
    #self.cu[3] = -1.0;
    #self.cu[4] = -1.0;
    self.rhs = _rhs   
 
    ## derivative operators 
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
    self.K = np.asfortranarray(
              self.Ka1b1c_a1b1c1.dot(self.Ka1bc_a1b1c.dot(self.Kabc_a1bc)))
    # promotion so derivs are in same basis
    self.K_a1bc1_a1b1c1 = np.asfortranarray(
                            kMat(_a+1, _b, _c+1, self.Ha1bc1, self.Ha1b1c1, _n-1, 1))
    self.K_ab1c1_a1b1c1 = np.asfortranarray(
                            kMat(_a, _b+1, _c+1, self.Hab1c1, self.Ha1b1c1, _n-1, 0))

    self.Dx = np.asfortranarray(
              self.K_a1bc1_a1b1c1.dot(
                dMat  ( self.a, self.b, self.c,
                        self.Habc, self.Ha1bc1, _n-1, 0 ) ) )
    self.Dy = np.asfortranarray(
              self.K_ab1c1_a1b1c1.dot(
                dMat  ( self.a, self.b, self.c,
                        self.Habc, self.Hab1c1, _n-1, 1 ) ) )
    self.rhscoeffs = np.asfortranarray(
                      self.K.dot( 
                        self.rhsCoeffs() ) )


    ## boundary evaluators
    self.nl = self.N
    self.nb = self.N
    self.nh = self.N
    self.xl = np.zeros((self.nl,))
    self.yl, _ = legendre(self.nl); self.yl = (self.yl + 1.0) / 2.0 
    self.xb, _ = legendre(self.nb); self.xb = (self.xb + 1.0) / 2.0 
    self.yb = np.zeros((self.nb,))
    self.xh, _ = legendre(self.nh); self.xh = (self.xh + 1.0) / 2.0 
    self.yh = 1.0 - self.xh
    self.polyl = jPoly(self.nl, _n, _a, _b, _c, 1)
    self.polyb = jPoly(self.nb, _n, _a, _b, _c, 1)
    self.polyh = jPoly(self.nh, _n, _a, _b, _c, 1)
    self.vl = computeV(self.polyl, self.N, self.xl, self.yl)
    self.vb = computeV(self.polyb, self.N, self.xb, self.yb)
    self.vh = computeV(self.polyh, self.N, self.xh, self.yh)

    ## post / path processing
    
    # quadrature nodes
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
    self.polypt = jPoly(1, _n, _a, _b, _c, 1)
    self.nthetas = _nthetas
    self.polyCirc = jPoly(_nthetas, _n, _a, _b, _c, _nthreads)
    self.thetas = np.linspace(0, 2*np.pi, _nthetas) 
    self.r = _r
    self.Xcirc = np.asfortranarray(_r * np.cos(self.thetas))
    self.Ycirc = np.asfortranarray(_r * np.sin(self.thetas))


  def solveH(self):
    libeikonal.eikonalSolveH(self.N,
                             self.Dx.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)), 
                             self.Dy.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)),
                             self.rhscoeffs.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)), 
                             self.cu.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)))
    return self.cu 

  def solveP(self, ul, ub, uh):
    libeikonal.eikonalSolveP(self.N,
                             self.nl,
                             self.nb,
                             self.nh, 
                             self.Dx.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)), 
                             self.Dy.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)),
                             self.rhscoeffs.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)), 
                             self.vl.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)),
                             ul.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)),
                             self.vb.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)), 
                             ub.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)),
                             self.vh.ctypes.data_as(
                               ctypes.POINTER(
                               ctypes.c_double)), 
                             uh.ctypes.data_as(
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
   
    Qrule = np.loadtxt("triquadLeg_10_17.txt"); 
    N = int(len(Qrule) / 3)
    X = Qrule[0:N]
    Y = Qrule[N:2*N]
    W = Qrule[2*N:3*N]
    cF, _ = computeCoeffs(self.rhs, self.n-1, self.a, self.b, self.c, X, Y, W, 1)
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
  
  def cuCoeffs(self):
   
    Qrule = np.loadtxt("triquadLeg_10_17.txt"); 
    N = int(len(Qrule) / 3)
    X = Qrule[0:N]
    Y = Qrule[N:2*N]
    W = Qrule[2*N:3*N]
    print(self.distToTri(X,Y))
    cu, _ = computeCoeffs(self.distToTri, self.n-1, self.a, self.b, self.c, X, Y, W, 1)
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
                                     ctypes.POINTER(ctypes.c_double)]
libeikonal.eikonalSolveH.restype = None                                    

libeikonal.eikonalSolveP.argtypes = [ctypes.c_uint,\
                                     ctypes.c_uint,\
                                     ctypes.c_uint,\
                                     ctypes.c_uint,\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double),\
                                     ctypes.POINTER(ctypes.c_double)]
libeikonal.eikonalSolveP.restype = None                                    
