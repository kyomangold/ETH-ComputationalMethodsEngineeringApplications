from numpy import *
import matplotlib
import matplotlib.pyplot as plt

error = loadtxt('error.txt')
dt = loadtxt('dt.txt')

# We ignore the first few datapoints,
# as we can see from the plot they are a bit wrong
rate = polyfit(log(dt[2:]), log(error[2:]), 1)[0]
print rate

plt.loglog(dt, error, '-o', label='$|u(T)-U_N|$')
plt.loglog(dt, error[-1]*dt**rate/(dt[-1]**rate), '--', label='$\Delta t^{%0.2f}$' % rate)
plt.xlabel('$\Delta t$')
plt.ylabel('$Error$')
#plt.xlim([dt[0], dt[1]])
plt.legend(loc="best")
plt.show()
