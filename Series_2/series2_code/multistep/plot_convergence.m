clc
close all
clear all

FolderName = 'build';

errors_euler= load([FolderName,'/errors.txt']);
resolutions = load([FolderName,'/resolutions.txt']);

loglog(resolutions, errors_euler);
xlabel('resolutions');
ylabel('errors');

if exist([FolderName,'/errors_ab2.txt'], 'file') == 2 % exists and is a file
    hold on
    errors_ab2 = load([FolderName,'/errors_ab2.txt']);
    loglog(resolutions, errors_ab2);
    legend('euler', 'ab2');
else
    legend('euler');
end
grid minor