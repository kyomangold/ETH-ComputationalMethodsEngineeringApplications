from numpy import *
from pylab import *

classes = ["S","E","I","R","D"]

t = loadtxt('time.txt')

for c in classes :
	y = loadtxt(c + ".txt")
	plot(t,y,label = c)

xlabel(r'$Time$')
ylabel(r'$Individuals$')

grid(True)
legend(loc = "best")

# savefig("../../figs/SEIR_fe_" + str(t.size-1) + ".png")
# savefig("../../figs/SEIR_dirk_" + str(t.size-1) + ".png")

show()

