clc
close all
clear all

FolderName = 'build';

y = load([FolderName,'/y.txt']);
t = load([FolderName,'/time.txt']);
plot(t,y, 'DisplayName', 'y');

legend TOGGLE