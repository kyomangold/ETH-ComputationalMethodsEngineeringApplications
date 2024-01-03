function M=sol_movie(sol_filename)
u=load(sol_filename);
t=load(strcat("time", extractAfter(sol_filename, 1)));
x = linspace(0,5,size(u,2));
Nf = size(u,1); % number of frames
j=1;

axis([0 5 0 2.1]);
xlabel('x','Fontsize',14)
ylabel('u(x)','Fontsize',14)
hold on

for i=1:Nf
    cla
    plot(x,u(i,:))
    title(strcat("T = ", num2str(t(i))));
    pause(0.001)
    M(j)=getframe(gcf);
    j=j+1;
end

