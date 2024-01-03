from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import sys

error = loadtxt('errors.txt')
numbers = loadtxt('numbers.txt')
times = loadtxt('walltimes.txt')

n_resolutions = 10

if max(error) == 0 or max(numbers) == 0:
	print("numbers or errors is identically zero!")
	sys.exit()

# We ignore the first few datapoints,
# as they may be underresolved
rate = polyfit(log(numbers[2:]), log(error[2:]), 1)[0]
print("Error scales with number of cells with order", rate)
plt.loglog(numbers, error, '-o', label='$|u(T)-U_N|$')

if (error.size < n_resolutions):
	plt.savefig("../../figs/error_vs_N.png")
	plt.loglog(numbers, error[-1]*numbers**rate/(numbers[-1]**rate), '--', label='$\Delta t^{%0.2f}$' % rate)
else :
	plt.savefig("../../figs/error_vs_N_satur.png")

plt.xlabel('$N$')
plt.ylabel('$Error$')
#plt.xlim([dt[0], dt[1]])
plt.legend(loc=0)
plt.grid(True)

plt.show()


if max(times) > 0:
	rate = polyfit(log(numbers[2:]), log(times[2:]), 1)[0]
	print("Runtime scales with number of cells with order", rate)
	plt.loglog(numbers, times, '-o', label='Walltimes')

	if (error.size < n_resolutions):
		plt.savefig("../../figs/walltime_vs_N.png")
		plt.loglog(numbers, times[-1]*numbers**rate/(numbers[-1]**rate), '--', label='$\Delta t^{%0.2f}$' % rate)
	else :
		plt.savefig("../../figs/walltime_vs_N_satur.png")

	plt.xlabel('$N$')
	plt.ylabel('$Time$')
	plt.legend(loc=0)
	plt.grid(True)

	plt.show()


	rate = polyfit(log(times[2:]), log(error[2:]), 1)[0]
	print("Error scales with runtime with order", rate)
	plt.loglog(times, error, '-o', label='$|u(T)-U_N|$')

	if (error.size < n_resolutions):
		plt.savefig("../../figs/error_vs_walltime.png")
		plt.loglog(times, error[-1]*times**rate/(times[-1]**rate), '--', label='$\Delta t^{%0.2f}$' % rate)
	else :
		plt.savefig("../../figs/error_vs_walltime_satur.png")

	plt.xlabel('$Time$')
	plt.ylabel('$Error$')
	plt.legend(loc=0)
	plt.grid(True)

	plt.show()
