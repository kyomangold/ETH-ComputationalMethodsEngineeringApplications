import numpy, pylab
u = numpy.loadtxt("u_wo_bc_fem.txt")
x = numpy.loadtxt("x_wo_bc_fem.txt")
pylab.plot(x, u, '-o', label='$u(x)$')
pylab.plot(x, numpy.sin(x), label='$\\sin(x)$')
pylab.legend()
pylab.show()

ubc = numpy.loadtxt("u_w_bc_fem.txt")
xbc = numpy.loadtxt("x_w_bc_fem.txt")
pylab.plot(xbc, ubc, '-o', label='$u(x)$')
pylab.plot(xbc, numpy.cos(xbc)+ 3./(4*numpy.pi)*xbc+3.0/4.0, label='$\\cos(x)+\\frac{3}{4\\pi}x + \\frac{3}{4}$')
pylab.legend()
pylab.show()
