from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation
import sys



class HeatDataAnimator(object):
    def __init__(self, u, method, bc, ts, a):
        self.u = u
        self.x = linspace(0,1,u.shape[1])

        self.fig = plt.figure()

        self.ax = plt.axes(xlim=(-0.1, 1.1), ylim=(amin(u)-0.1, amax(u)+0.1))
        self.ax.grid(True)
        self.ax.set_title("{method} with ${g}_L, {g}_R$, $\\Delta t={ts}$, $a={a}$".\
                       format(method=method, g=bc, ts=ts, a=a))
        self.line, = self.ax.plot([], [], lw=2)
        self.line.set_data([], [])

        self.time_title = self.ax.text(0.2,0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=self.ax.transAxes, ha="center")
        self.ts = ts

    def __call__(self, i):
        self.time_title.set_text("Current time: {:.3f}".format(i*self.ts))
        self.line.set_data(self.x, self.u[i,:])
        return self.line, self.time_title


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage:\n\tpython {} <method>\n".format(sys.argv[0]))
        print("Where method is either forward_euler or crank_nicolson")
        exit(1)

    method = sys.argv[1]
    boundary_conditions = ['g^1', 'g^2', 'g^3']

    for i, bc in enumerate(boundary_conditions):

        filename = '{method}_boundaries_u_g{i}.txt'.format(
            i=i, method=method)

        u = loadtxt(filename)


        animator = HeatDataAnimator(u, method, bc, 0.5/64**2, 'a_1')

        anim = matplotlib.animation.FuncAnimation(animator.fig,
                                                  animator,
                                                  frames=u.shape[0],
                                                  interval=5,
                                                  blit=True)

        plt.show()
