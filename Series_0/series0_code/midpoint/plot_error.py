import numpy
import pylab
errors = numpy.loadtxt('series0_1_d_errors.txt')
numberOfIntervals = numpy.loadtxt('series0_1_d_numbers.txt')

pylab.loglog(numberOfIntervals, errors, '-*')
pylab.title('Error plot for ex. 0.1d)')
pylab.xlabel('n')
pylab.ylabel('|I(f)-I_n(f)|')
pylab.grid('on')
pylab.show()
