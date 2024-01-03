clc %to clear the command
close all % to close eventual previous windows
clear all % to clear the workspace

FolderName = 'build/';

S = load([FolderName,'S.txt']);
E = load([FolderName,'E.txt']);
I = load([FolderName,'I.txt']);
R = load([FolderName,'R.txt']);
D = load([FolderName,'D.txt']);
t = load([FolderName,'time.txt']);

plot(t, S, 'b', t, E, 'r', t, I, 'g', t, R, 'k', t, D, 'm')
xlabel("Time");
ylabel("Individuals");

legend('S', 'E','I', 'R', 'D')