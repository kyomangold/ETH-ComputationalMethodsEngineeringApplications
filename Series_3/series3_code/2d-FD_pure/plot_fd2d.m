u = load("u_fd.txt");
N = sqrt(length(u));
disp("N=" + N)
u = reshape(u,[N,N])';

clear title xlabel ylabel zlabel
[x,y] = meshgrid(linspace(0,1,N), linspace(0,1,N));
surf(x,y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title("Approximation to the solution")
disp('Press any key to continue...')
pause
cla

% clear title xlabel ylabel zlabel
% surf(x,y, sin(2*pi*x).*cos(2*pi*y))
% xlabel('x');
% ylabel('y');
% zlabel('u');
% disp('Press any key to continue...')
% pause
% cla

clear title xlabel ylabel zlabel
contourf(x,y,abs(u - sin(2*pi*x).*sin(2*pi*y)))
colorbar()
title("Approximation error")
disp('Press any key to continue...')
pause
cla


clear title xlabel ylabel zlabel
n = floor(N/2);
plot(x(n,:), u(n,:));
hold on
plot(x(n,:), sin(2*pi*x(n,:)).*sin(2*pi*y(n,:)));
legend('u', 'exact');
xlabel = 'x';
ylabel = 'u';
title("Cut at y = " + n/(N+1) );
