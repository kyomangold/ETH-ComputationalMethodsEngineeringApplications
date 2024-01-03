from numpy import loadtxt, exp, sqrt, sin, cos
import pylab as plt

time_expl=loadtxt("time_expl.txt");
pos_expl=loadtxt("position_expl.txt");
v_expl=loadtxt("velocity_expl.txt");
en_expl=loadtxt("energy_expl.txt");
time_impl=loadtxt("time_impl.txt");
pos_impl=loadtxt("position_impl.txt");
v_impl=loadtxt("velocity_impl.txt");
en_impl=loadtxt("energy_impl.txt");

zeta=0.2;
freq = sqrt(1-zeta**2);
exact_pos = exp(-zeta*time_impl) * (zeta/freq*sin(freq*time_impl)+cos(freq*time_impl));

plt.plot(time_expl,pos_expl,'b', label='explicit')
plt.plot(time_impl,pos_impl,'r', label='implicit')
plt.plot(time_impl,exact_pos,'g', label='exact')
plt.title('Position')
plt.legend()
plt.xlabel('t')
plt.ylabel('x(t)')
plt.show()

plt.plot(time_expl,v_expl, label='explicit')
plt.plot(time_impl,v_impl,'r', label='implicit')
plt.title('Velocity')
plt.legend()
plt.xlabel('t')
plt.ylabel('v(t)')
plt.show()

plt.plot(time_expl,en_expl, label='explicit')
plt.plot(time_impl,en_impl,'r', label='implicit')
plt.title('Kinetic energy')
plt.legend()
plt.xlabel('t')
plt.ylabel('E(t)')
plt.show()
