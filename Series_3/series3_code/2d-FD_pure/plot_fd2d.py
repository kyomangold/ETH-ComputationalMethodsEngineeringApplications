import mpl_toolkits.mplot3d
from numpy import *
from pylab import *

u = loadtxt("u_fd.txt")
N = int(sqrt(u.shape[0]))
print("N=%d" % N)
u.shape = (N,N)

fig = figure()
ax = fig.gca(projection="3d")
x,y = mgrid[0:1:N*1j, 0:1:N*1j]
ax.plot_surface(x,y,u)
xlabel('x')
ylabel('y')
#zlabel('u')
title('Approximation to the solution')
show()

pcolormesh(x,y,abs(u - sin(2*pi*x)*sin(2*pi*y)))
colorbar()
xlabel('x')
ylabel('y')
title('Approximation error')
show()

n = int(N/2)
plot(x[:, n], u[:, n], label='u')
plot(x[:,n], sin(2*pi*x[:,n])*sin(2*pi*y[:,n]), label='exact')
xlabel('x')
ylabel('u')
title('Cut at y={}'.format(round(n*1./(N+1), 5)))
legend()
show()
