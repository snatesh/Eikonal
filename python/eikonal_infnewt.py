from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
# Create mesh and function space

#domain = Polygon([Point(0,0), Point(1,0), Point(0,1)])
#mesh = generate_mesh(domain, 64)
#V = FunctionSpace(mesh,'CG',1)
#class DirichletBoundary(SubDomain):
#  def inside(self, x, on_boundary):
#    return near(x[0],0) or near(x[1],0) or \
#           near(x[0]+x[1],1)
domain = Circle(Point(0,0),1)
mesh = generate_mesh(domain,64)

#Defining boundary condition
V = FunctionSpace(mesh,'CG',1)

#class DirichletBoundary(SubDomain):
#    def inside(self, x, on_boundary):
#        return near(x[0], 0) or near(x[0], 1.0) or \
#               (near(x[1],0) or near(x[1],1.0))

class DirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
	  if (x[0]*x[0] + x[1]*x[1] > 0.99):
	  	return True
	  else:
	  	return False

u0 = Constant(0.0)
u0_boundary = DirichletBoundary()
bc = DirichletBC(V, u0, u0_boundary, "pointwise")

# Define functions
u = Function(V)  # Solution function
du = Function(V)  # Newton step (correction)
v = TestFunction(V)
w = TrialFunction(V)

# Define parameters
xi = Constant(0.01)  # Parameter in the functional
f = Constant(1.0)

# Define the energy functional
E = inner(grad(u), grad(u))**0.5 * dx + (xi/2) * inner(grad(u), grad(u)) * dx - f * u * dx

# Compute the first variation (Residual F = dE/du)
F = derivative(E, u, v)

# Compute the second variation (Hessian H = d²E/du²)
H = derivative(F, u, w)

# Assemble the system
#u_0 = Expression('x[0]*x[1]*(1-x[0]-x[1])',degree=1)
#u.assign(interpolate(u_0, V))
F1 = inner(grad(w), grad(v))*dx - f*v*dx
a1, L1 = lhs(F1), rhs(F1)
solve(a1==L1, u, bc)

rtol = 1e-8
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
  myit = solve(J, du.vector(), b, "cg", "petsc_amg") 
  lin_it = lin_it + myit
  # Update solution u <- u + du
  u.vector().axpy(-1., du.vector())
  pi = assemble(E)
  print("{0:3d} {1:3d} {2:15e} {3:15e} {4:15e}".format(
        i, myit, pi, -b.inner(du.vector()), b.norm("l2")))

XY = V.tabulate_dof_coordinates()
X = XY[:,0]
Y = XY[:,1]

plt.figure()
plot(u)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(X,Y,u.vector().vec())

### compute initial guess for newton iter
F1 = inner(grad(w), grad(v))*dx - f*v*dx
a1, L1 = lhs(F1), rhs(F1)
solve(a1==L1, u, bc)

F = sqrt(inner(grad(u), grad(u)))*v*dx -\
    f*v*dx + xi*inner(grad(u),grad(v))*dx
solve(F==0, u, bc)
plt.figure()
plot(u)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(X,Y,u.vector().vec())

plt.show()


