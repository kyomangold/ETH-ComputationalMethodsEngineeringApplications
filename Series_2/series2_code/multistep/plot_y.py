from numpy import *
from pylab import *

y = loadtxt("y.txt")
t = loadtxt('time.txt')
plot(t, y, label="y")
legend(loc = "best")
show()

