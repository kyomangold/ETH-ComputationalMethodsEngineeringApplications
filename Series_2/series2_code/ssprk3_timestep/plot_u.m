clc
close all
clear all

FolderName = 'build';
u = load([FolderName,'/solution.txt']);
t = load([FolderName,'/time.txt']);

plot(t, u);
xlabel('t');
ylabel('u(t)');

