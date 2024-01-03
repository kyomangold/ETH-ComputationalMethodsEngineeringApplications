u1 = load('u1.txt');
u2 = load('u2.txt');
u3 = load('u3.txt');
u4 = load('u4.txt');
t = load('time.txt');

plot(t, u1, t, u2, t, u3, t, u4)
legend('u_1', 'u_2', 'u_3', 'u_4')
figure()
plot(t, u1+u2+u3+u4)
legend('u_1+u_2+u_3+u_4')
