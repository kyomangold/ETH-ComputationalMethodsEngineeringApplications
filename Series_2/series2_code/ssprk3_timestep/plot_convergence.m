clc
close all
clear all

FolderName = 'build';
error = load([FolderName,'/error.txt']);
dt = load([FolderName,'/dt.txt']);

N = length(error);
rate = polyfit(log(dt(2:N)), log(error(2:N)), 1);
rate = rate(1)

loglog(dt, error, '-o')
hold on
loglog(dt, error(N-1)*dt.^rate/(dt(N-1).^rate), '--')
xlabel('\Delta t')
ylabel('Error')
axis([dt(N) dt(1) -inf inf]);
label = strcat('\Delta t^{', num2str(rate,3), '}');
legend('|u(T)-U_N|', label, 'Location', 'northwest')
grid on
