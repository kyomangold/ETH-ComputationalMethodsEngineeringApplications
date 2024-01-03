import numpy
import matplotlib
import matplotlib.pyplot as plt
errors = numpy.loadtxt('series0_1_d_errors.txt')
numberOfIntervals = numpy.loadtxt('series0_1_d_numbers.txt')

plt.loglog(numberOfIntervals, errors, '-*', label='$|I_n(f)-I(f)|$')
plt.title('Error plot for the midpoint rule')
plt.xlabel('$n$')
plt.ylabel('$|I(f)-I_n(f)|$')
plt.grid('on')

# Now compare against the rate n^(-2)
plt.hold('on')
plt.loglog(numberOfIntervals, numberOfIntervals**(-2), '--',label='$\mathcal{O}(n^{-2})$')
plt.legend()
plt.show()
