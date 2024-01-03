u = load("u_wo_bc_fem.txt");
x = load("x_wo_bc_fem.txt");
plot(x, u, '-o');
hold on
plot(x, sin(x))
leg = legend('$u(x)$', '$\\sin(x)$');
set(leg, 'Interpreter', 'latex');
set(leg, 'FontSize', 12);
xlabel("x");
ylabel("u");
title("Without BC");
disp('Press any key to continue...')
pause
cla

ubc = load("u_w_bc_fem.txt");
xbc = load("x_w_bc_fem.txt");
plot(xbc, ubc, '-o');
hold on
plot(xbc, cos(xbc)+ 3./(4*pi)*xbc+3.0/4.0);
leg = legend('$u(x)$', '$\\cos(x)+\frac{3}{4\pi}x + \frac{3}{4}$');
set(leg, 'Interpreter', 'latex');
set(leg, 'FontSize', 12);
xlabel("x");
ylabel("u");
title("Without BC");
