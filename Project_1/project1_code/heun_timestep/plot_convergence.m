clc
close all
clear all

FolderName = 'build/';
error = load([FolderName,'error.txt']);
dt = load([FolderName,'dt.txt']);

rates = polyfit(log(dt(2:end)),log(error(2:end)),1);
rate = rates(1);

P = loglog(dt, error, '-o');
hold on
Pr = loglog(dt,(error(end)/dt(end)^rate)*(dt).^rate,'--');

Ps = [P, Pr];
Ps_labs = {'$\vert u(T)-U_N \vert$',['$',strcat('\Delta t^{',num2str(rate)),'}$']};

legend(Ps,Ps_labs,'interpreter','latex','Location', 'NorthWest','FontSize',13);

xlabel('\Delta t');
ylabel('Error');

