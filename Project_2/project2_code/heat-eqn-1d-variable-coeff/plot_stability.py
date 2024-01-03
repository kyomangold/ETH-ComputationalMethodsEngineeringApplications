import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage:\n\tpython {} <method name>\n".format(sys.argv[0]))
        print("where <method name> is either forward_euler or crank_nicolson")
        exit(1)

    method = sys.argv[1]

    coefficients = ['a_1', 'a_2', 'a_3']

    N = 127
    dxSquared = 1.0 / (N+1)**2

    timesteps = 0.5*dxSquared * np.array([128., 8., 1.])
    fig, axes = plt.subplots(len(coefficients), len(timesteps))

    for i, coefficient in enumerate(coefficients):
        for j, timestep in enumerate(timesteps):
            filename = '{method}_stability_u_a{i}_dt{j}.txt'.format(
                i=i, j=j, method=method)

            u = np.loadtxt(filename)
            x = np.linspace(0, 1, u.shape[1])

            axes[i,j].plot(x, u[-1,:])
            axes[i,j].set_title('$a_{i}$, $\\Delta t_{j}$'.format(i=i+1,
                                                                  j=j+1))
    fig.tight_layout()
    plt.savefig('stability_{}.png'.format(method))
    plt.show()
