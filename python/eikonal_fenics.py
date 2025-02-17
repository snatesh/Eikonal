from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
from convToSurf3D import surfTo3D

# Create mesh and define function space
#nx = 20; ny = 20
#mesh=RectangleMesh(Point(0, 0), Point(1, 1), nx, ny)
#domain = Rectangle(Point(0,0),Point(1,1))

#domain = Circle(Point(0,0),1)
#mesh = generate_mesh(domain,64)

domain = Polygon([Point(0,0), Point(1,0), Point(0,1)])
mesh = generate_mesh(domain, 32)

#Defining boundary condition
V = FunctionSpace(mesh,'CG',3)
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
u_hat = TestFunction(V)
u_tilde = TrialFunction(V)
xi = Constant(0.01)

f = Constant(1.0)
n = FacetNormal(mesh)

# Functional we seek to minimize
#Pi = (sqrt(inner(nabla_grad(u), nabla_grad(u)) + xi**2*div(nabla_grad(u))**2) - f)*dx
#Pi = (sqrt(inner(nabla_grad(u), nabla_grad(u))) - xi*div(nabla_grad(u)) - f)*dx
#Pi = (sqrt(inner(nabla_grad(u),nabla_grad(u)))*u - f*u + xi*inner(nabla_grad(u),nabla_grad(u)))*dx
#Pi = f*sqrt(inner(nabla_grad(u),nabla_grad(u)))*dx
#Pi = (sqrt(inner(nabla_grad(u),nabla_grad(u))) + (xi/2.0)*inner(nabla_grad(u),nabla_grad(u)) - f*u)*dx
#Pi = (Constant(1./3.)*sqrt(inner(nabla_grad(u),nabla_grad(u)))**3 - \
#        0.5*f*u + (xi/2.0)*div(nabla_grad(u))**2)*dx
#Pi = (Constant(0.5)*f*inner(grad(u),grad(u)) + (xi/0.5)*div(grad(u))**2 - Constant(0.5))*dx
#Pi = np.power((inner(grad(u),grad(u)) - (1.0/f**2) + xi*div(grad(u)))**2,0.25)*dx
#Pi = (sqrt(inner(grad(u),grad(u))) - f*u - xi*div(grad(u)))**2*dx
Pi = (sqrt(inner(grad(u),grad(u))) + xi*inner(grad(u),grad(u)) - f*u)*dx
# first variation of functional
#dPi = ((inner(nabla_grad(u),nabla_grad(u_hat)) + xi**2*div(nabla_grad(u_hat))*div(nabla_grad(u))) / \
#			 sqrt(inner(nabla_grad(u),nabla_grad(u)) + xi**2*div(nabla_grad(u))**2))*dx
#dPi = (inner(nabla_grad(u),nabla_grad(u_hat)) / \
#			 sqrt(inner(nabla_grad(u),nabla_grad(u))) - xi*div(nabla_grad(u_hat)) )*dx
#dPi = u*(inner(nabla_grad(u),nabla_grad(u_hat)) / \
#         sqrt(inner(nabla_grad(u),nabla_grad(u))))*dx + \
#      u_hat*(sqrt(inner(nabla_grad(u),nabla_grad(u))))*dx - \
#      f*u_hat*dx + xi*(2*inner(nabla_grad(u),nabla_grad(u_hat)))*dx

#dPi = (f*inner(nabla_grad(u),nabla_grad(u_hat)) / \
#       sqrt(inner(nabla_grad(u),nabla_grad(u))))*dx
#dPi = (inner(nabla_grad(u),nabla_grad(u_hat))/(sqrt(inner(nabla_grad(u),nabla_grad(u)))) +\
#      xi*inner(nabla_grad(u),nabla_grad(u_hat)) - f*u_hat)*dx
dPi = derivative(Pi, u, u_hat)
# verify expression for first variation
u0 = interpolate(Expression("x[0]*(x[0]-1)*x[1]*(x[1]-1)",degree=1), V)
n_eps = 32
eps = 1e-1*np.power(2., -np.arange(n_eps))
err_grad = np.zeros(n_eps)

u.assign(u0)
pi0 = assemble(Pi)
grad0 = assemble(dPi)

direction = Function(V)
direction.vector().set_local(np.random.randn(V.dim()))
direction_grad0 = grad0.inner(direction.vector())

for i in range(n_eps):
    u.assign(u0)
    u.vector().axpy(eps[i], direction.vector()) #u = u + eps[i]*direction
    piplus = assemble(Pi) 
    err_grad[i] = abs( (piplus - pi0)/eps[i] - direction_grad0 )
print(err_grad)
plt.figure()
plt.loglog(eps, err_grad, "-ob")
plt.loglog(eps, (.5*err_grad[0]/eps[0])*eps, "-.k")
plt.title("Finite difference check of the first variation (gradient)")
plt.xlabel("eps")
plt.ylabel("Error grad")
plt.legend(['Error Grad', 'First Order'], loc='upper left')
plt.show()
# second variation of functional
#ddPi = ((inner(nabla_grad(u_tilde),nabla_grad(u_hat))+\
#         xi**2*div(nabla_grad(u_tilde))*div(nabla_grad(u_hat))) / \
#        sqrt(inner(nabla_grad(u),nabla_grad(u))+xi**2*div(nabla_grad(u))**2))*dx - \
#       ((inner(nabla_grad(u),nabla_grad(u_hat))+xi**2*div(nabla_grad(u))*div(nabla_grad(u_hat)))* \
#        (inner(nabla_grad(u),nabla_grad(u_tilde))+xi**2*div(nabla_grad(u))*div(nabla_grad(u_tilde))) / \
#        np.power(inner(nabla_grad(u),nabla_grad(u))+xi**2*div(nabla_grad(u))**2,1.5))*dx
#
#ddPi = (inner(nabla_grad(u_tilde),nabla_grad(u_hat))/sqrt(inner(nabla_grad(u),nabla_grad(u))))*dx - \
#       (inner(nabla_grad(u),nabla_grad(u_tilde))*inner(nabla_grad(u),nabla_grad(u_hat)) / \
#       np.power(inner(nabla_grad(u),nabla_grad(u)),1.5))*dx + \
#       xi*inner(nabla_grad(u_tilde),nabla_grad(u_hat))*dx
ddPi = derivative(dPi, u, u_tilde)
# verify expression for second variation
u.assign(u0)
H_0 = assemble(ddPi)
err_H = np.zeros(n_eps)
for i in range(n_eps):
    u.assign(u0)
    u.vector().axpy(eps[i], direction.vector())
    grad_plus = assemble(dPi)
    diff_grad = (grad_plus - grad0)
    diff_grad *= 1/eps[i]
    H_0dir = H_0 * direction.vector()
    err_H[i] = (diff_grad - H_0dir).norm("l2")

plt.figure()
plt.loglog(eps, err_H, "-ob")
plt.loglog(eps, (.5*err_H[0]/eps[0])*eps, "-.k")
plt.title("Finite difference check of the second variation (Hessian)")
plt.xlabel("eps")
plt.ylabel("Error Hessian")
plt.legend(['Error Hessian', 'First Order'], loc="upper left")
plt.show()

# Starting guess satisfies Dirichlet bdry condition


# compute initial guess for newton iter
#F1 = inner(grad(u_tilde), grad(u_hat))*dx - f*u_hat*dx
#a1, L1 = lhs(F1), rhs(F1)
#solve(a1==L1, u, bc)
# solve grad==0 given initial guess

u_0 = Expression('x[0]*x[1]*(1-x[0]-x[1])',degree=1)
u.assign(interpolate(u_0, V))
parameters={"newton_solver": {"relative_tolerance": 1e-10,\
                               "report": True, \
                               "maximum_iterations": 100}}
solve(dPi == 0, u, bc, solver_parameters=parameters)
print("Built-in FEniCS non linear solver.")
print("Norm of the gradient at converge", assemble(dPi).norm("l2"))
print("Value of the energy functional at convergence", assemble(Pi))
XY = V.tabulate_dof_coordinates()
X = XY[:,0]
Y = XY[:,1]

plt.figure()
plot(u)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(X,Y,u.vector().vec())
plt.show()
#
#
#
##
### compute initial guess for newton iter
##F1 = inner(grad(u_tilde), grad(u_hat))*dx - f*u_hat*dx
##a1, L1 = lhs(F1), rhs(F1)
##solve(a1==L1, u, bc)
##
##F = sqrt(inner(grad(u), grad(u)))*u_hat*dx -\
##    f*u_hat*dx + xi*inner(grad(u),grad(u_hat))*dx
##solve(F==0, u, bc)
##
##
##XY = V.tabulate_dof_coordinates()
##X = XY[:,0]
##Y = XY[:,1]
##
##plt.figure()
##plot(u)
##fig = plt.figure()
##ax = fig.add_subplot(projection='3d')
##ax.scatter(X,Y,u.vector().vec())
##plt.show()
