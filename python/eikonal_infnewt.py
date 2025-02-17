from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
# Create mesh and function space

domain = Polygon([Point(0,0), Point(1,0), Point(0,1)])
mesh = generate_mesh(domain, 64)
V = FunctionSpace(mesh,'CG',1)
XY = V.tabulate_dof_coordinates()
X = XY[:,0]
Y = XY[:,1]
class DirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0],0) or near(x[1],0) or \
           near(x[0]+x[1],1)
#domain = Circle(Point(0,0),1)
#mesh = generate_mesh(domain,64)

#class DirichletBoundary(SubDomain):
#    def inside(self, x, on_boundary):
#        return near(x[0], 0) or near(x[0], 1.0) or \
#               (near(x[1],0) or near(x[1],1.0))

#class DirichletBoundary(SubDomain):
#  def inside(self, x, on_boundary):
#	  if (x[0]*x[0] + x[1]*x[1] > 0.99):
#	  	return True
#	  else:
#	  	return False

u0 = Constant(0.0)
u0_boundary = DirichletBoundary()
bc = DirichletBC(V, u0, u0_boundary, "pointwise")

# Define functions
u = Function(V)  # Solution function
du = Function(V)  # Newton step (correction)
v = TestFunction(V)
w = TrialFunction(V)

# Define parameters
xi = Constant(0.001)  # Parameter in the functional
f = Constant(2.0)
#f = Expression('x[0]*x[1]*(1-x[0]-x[1])',degree=1)

# proxy for energy
E = ((sqrt(inner(grad(u),grad(u))) - xi*div(grad(u)) - (1./f)))*dx

# Compute the first variation (Residual F = dE/du)
F = (sqrt(inner(grad(u),grad(u)))*v - (1./f)*v + xi*inner(grad(u),grad(v)))*dx
#F = inner(grad(u),grad(u))*v*dx - v/f**2*dx - xi**2*div(grad(u))*div(grad(v))*dx



# Compute the second variation (Hessian H = d²E/du²)
H = derivative(F, u, w)
#H = (inner(grad(u),grad(w))*v/(sqrt(inner(grad(u),grad(u)))))*dx + xi*inner(grad(w),grad(v))*dx

# Assemble the system
#u_0 = Expression('x[0]*x[1]*(1-x[0]-x[1])',degree=1)
#u.assign(interpolate(u_0, V))

F1 = inner(grad(w), grad(v))*dx - (1./f)*v*dx
a1, L1 = lhs(F1), rhs(F1)
solve(a1==L1, u, bc)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(X,Y,u.vector().vec())

rtol = 1e-10
max_iter = 100

pi0 = assemble(E)
g0 = assemble(F)
tol = g0.norm("l2")*rtol
du = Function(V)

lin_it = 0
print("{0:3} {1:3} {2:15} {3:15} {4:15}".format(
      "It", "cg_it", "Energy", "(g,du)", "||g||l2"))

for i in range(max_iter):
  [J, b] = assemble_system(H, F, bc)
  if b.norm("l2") < tol:
    print("\nConverged in ", i, "Newton iterations and ", lin_it, "linear iterations.")
    break
  myit = solve(J, du.vector(), b) 
  lin_it = lin_it + myit
  # Update solution u <- u + du
  u.vector().axpy(-1., du.vector())
  pi = assemble(E)
  print("{0:3d} {1:3d} {2:15e} {3:15e} {4:15e}".format(
        i, myit, pi, -b.inner(du.vector()), b.norm("l2")))



plt.figure()
plot(u)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(X,Y,u.vector().vec())
plt.show()
