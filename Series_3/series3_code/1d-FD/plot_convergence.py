from numpy import *
from pylab import *
errors = loadtxt("errors_fd.txt")
resolutions = loadtxt("resolutions_fd.txt")
rate = polyfit(log(resolutions), log(errors), 1)[0]
print("Approximate order: %1.2f" %
      round(rate, 2))

loglog(resolutions, errors, '-o', label="$||u_N-u||_\\infty$")
loglog(resolutions, errors[-1]*resolutions**rate/(resolutions[-1]**rate), '--', label='$N^{%0.2f}$' % rate)
legend()
ylabel("errors")
xlabel("$N$")
grid("on")
show()

