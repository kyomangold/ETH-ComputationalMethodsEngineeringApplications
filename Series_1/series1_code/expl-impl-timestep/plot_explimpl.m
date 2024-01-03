time_expl=load('time_expl.txt','-ascii');
pos_expl=load('position_expl.txt','-ascii');
v_expl=load('velocity_expl.txt','-ascii');
en_expl=load('energy_expl.txt','-ascii');
time_impl=load('time_impl.txt','-ascii');
pos_impl=load('position_impl.txt','-ascii');
v_impl=load('velocity_impl.txt','-ascii');
en_impl=load('energy_impl.txt','-ascii');

zeta=0.2;
freq = sqrt(1-zeta^2);
exact_pos = exp(-zeta.*time_impl).*(zeta/freq*sin(freq*time_impl)+cos(freq*time_impl));

figure
plot(time_expl,pos_expl)
hold on
plot(time_impl,pos_impl,'r')
plot(time_impl,exact_pos,'g')
title('Position')
legend('explicit','implicit','exact')
xlabel('t')
ylabel('x(t)')

figure
plot(time_expl,v_expl)
hold on
plot(time_impl,v_impl,'r')
title('Velocity')
legend('explicit','implicit')
xlabel('t')
ylabel('v(t)')

figure
plot(time_expl,en_expl)
hold on
plot(time_impl,en_impl,'r')
title('Kinetic energy')
legend('explicit','implicit')
xlabel('t')
ylabel('E(t)')