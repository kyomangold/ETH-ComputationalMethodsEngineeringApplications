from numpy import *
import os
import matplotlib
import matplotlib.pyplot as plt

for dt in [2.5, 2, 1.5, 1, 0.5]:
    os.system("./heun %f" % dt)
    t = loadtxt('time.txt')
    u = loadtxt('solution.txt')
    plt.plot(t, u, '-*', label='$u(t)$')
    plt.xlabel('$t$')
    plt.ylabel('$u(t)$')
    plt.legend()
    plt.savefig('solution_%d.png' % int(10*dt))
    plt.close('all')
