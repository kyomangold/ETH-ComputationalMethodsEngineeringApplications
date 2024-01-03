import numpy
import pylab


u1 = numpy.loadtxt("u1.txt")
u2 = numpy.loadtxt("u2.txt")
u3 = numpy.loadtxt("u3.txt")
u4 = numpy.loadtxt("u4.txt")
t = numpy.loadtxt("time.txt")

pylab.plot(t, u1, label="u1")
pylab.plot(t, u2, label="u2")
pylab.plot(t, u3, label="u3")
pylab.plot(t, u4, label="u4")
pylab.legend()
pylab.show()

pylab.plot(t, u1+u2+u3+u4, label="u1+u2+u3+u4")
pylab.ylim([0,1.5])
pylab.legend()
pylab.show()
