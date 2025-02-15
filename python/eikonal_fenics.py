from dolfin import *
from mshr import *
import matplotlib.pyplot as plt


# Create mesh and define function space
#nx = 20; ny = 20
#mesh=RectangleMesh(Point(0, 0), Point(1, 1), nx, ny)
#domain = Rectangle(Point(0,0),Point(1,1))

#domain = Circle(Point(0,0),1)
#mesh = generate_mesh(domain,64)

domain = Polygon([Point(0,0), Point(1,0), Point(0,1)])
mesh = generate_mesh(domain, 64)

#Defining boundary condition
V = FunctionSpace(mesh,'CG',2)
u0 = Constant(0.0)

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

class DirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0],0) or near(x[1],0) or \
           near(x[0]+x[1],1)

u0_boundary = DirichletBoundary()
bc = DirichletBC(V, u0, u0_boundary, "pointwise")

u = Function(V)
du = TrialFunction(V)
v = TestFunction(V)
eps = Constant(0.001)

f = Constant(1.0)

# compute initial guess for newton iter
F1 = inner(grad(du), grad(v))*dx - f*v*dx
a1, L1 = lhs(F1), rhs(F1)
solve(a1==L1, u, bc)

F = sqrt(inner(grad(u), grad(u)))*v*dx - f*v*dx + eps*inner(grad(u),grad(v))*dx
solve(F==0, u, bc)


XY = V.tabulate_dof_coordinates()
X = XY[:,0]
Y = XY[:,1]

plt.figure()
plot(u)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(X,Y,u.vector().vec())
plt.show()
