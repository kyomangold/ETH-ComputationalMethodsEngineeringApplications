"""
Plots the linear regression and adds pretty cats to it all
"""

from matplotlib import pyplot as plt
import numpy 
import matplotlib
plt.xkcd()
Y=[ 6.5 , 7.2 , 7.3 , 7.6 , 7.7 , 7.9 , 7.9 , 7.9 , 7.9 , 8.0 , 8.3, 8.5 , 8.6 , 8.8 , 8.8]

X =[2 , 2.2 , 2.4 , 2.2 , 2.6 , 2.2 , 2.4 , 2.4 , 2.5 , 2.7 , 2.6 , 2.2 , 2.5 , 2.5 , 2.5 ]


plt.scatter(X, Y, 900*numpy.ones_like(X),c="g", alpha=0.5, marker=ur'$\U0001F431$', label="Observed data")


x = numpy.linspace(min(X), max(X),100)
regression = numpy.polyfit(X,Y, 1)
plt.plot(x, regression[0]*x+regression[1], '--', label="Regression")
plt.xlabel("Bodyweight (kg)")
plt.ylabel("Mass of heart (g)")

plt.legend(loc=0)
plt.show()


  
