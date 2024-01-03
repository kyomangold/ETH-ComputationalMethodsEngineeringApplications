from numpy import *
from pylab import *
from os.path import isfile

errors_euler= loadtxt("errors.txt")
resolutions = loadtxt("resolutions.txt")

loglog(resolutions, errors_euler, label="FE")
xlabel("resolutions")
ylabel("errors")
grid()

if isfile("errors_ab2.txt"):
	errors_ab2 = loadtxt("errors_ab2.txt")
	loglog(resolutions, errors_ab2, label="AB2")

legend()
show()

