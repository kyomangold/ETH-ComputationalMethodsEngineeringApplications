clear all
u = load("u_fd.txt");

u = [0 u 0]; % u contains interior points; add homogeneous BCs

x = linspace(0,1,length(u));
plot(x, u, '-o');
hold on
plot(x, (2*pi)^(-2) * sin(2*pi*x));
xlabel('x');
leg = legend('u', '$\frac{1}{(2\pi)^2}\sin(x)$');
set(leg, 'Interpreter', 'latex');
set(leg, 'FontSize', 12);
