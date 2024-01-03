clc
close all
clear all

FolderName = 'build/';
error = load([FolderName,'error.txt']);
resolutions = load([FolderName,'resolutions.txt']);

rates = polyfit(log(resolutions(2:end)),log(error(2:end)),1);
rate = rates(1);

P = loglog(resolutions, error, '-o');
hold on
Pr = loglog(resolutions,(error(end)/resolutions(end)^rate)*(resolutions).^rate,'--');

Ps = [P, Pr];
Ps_labs = {'$\vert u(T)-U_N \vert$',['$',strcat('\Delta t^{',num2str(rate)),'}$']};

legend(Ps,Ps_labs,'interpreter','latex','Location', 'NorthWest','FontSize',13);

xlabel('N');
ylabel('Error');