
from numpy import *
from pylab import *
import sys

basename = sys.argv[1]
errors = loadtxt("%s_errorsH1.txt" % basename)
numberOfElements = loadtxt("%s_resolutions.txt" % basename)
print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)

loglog(numberOfElements, errors, "-o")

xlabel("Number of free vertices")
ylabel("$||u_{FEM}-u||_{L^2}$")
grid("on")
show()

print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)
