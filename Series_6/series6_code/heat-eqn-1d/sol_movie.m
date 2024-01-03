function M=sol_movie(sol_filename)
u=load(sol_filename);
x = linspace(0,1,size(u,2));
Nf = size(u,1); % number of frames
j=1;

axis([0 1 min([u(1,:) -1]) max([1 u(1,:)])])
xlabel('x','Fontsize',14)
ylabel('u(x)','Fontsize',14)
hold on

for i=1:Nf
    cla
    plot(x,u(i,:))
    pause(0.001)
    M(j)=getframe(gcf);
    j=j+1;
end

