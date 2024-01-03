from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation
import sys



class HeatDataAnimator(object):
    def __init__(self, u):
        self.u = u
        self.x = linspace(0,1,u.shape[1])

        self.fig = plt.figure()
        self.ax = plt.axes(xlim=(0, 1), ylim=(0, 2))
        self.line, = self.ax.plot([], [], lw=2)
        self.line.set_data([], [])


    def __call__(self, i):
        self.line.set_data(self.x, self.u[i,:])
        return self.line,


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage:\n\tpython %s <u_*method*.txt>\n")
        exit(1)

    u = loadtxt(sys.argv[1])


    
    animator = HeatDataAnimator(u)

    anim = matplotlib.animation.FuncAnimation(animator.fig, animator,
                               frames=200, interval=20, blit=True)

    plt.show()

