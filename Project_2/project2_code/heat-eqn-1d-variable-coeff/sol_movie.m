function M=sol_movie(method_name)
    ulist = {};
    for i = 0:2
        sol_filename = sprintf('%s_boundaries_u_g%d.txt', method_name, i);
        u=load(sol_filename);
        ulist{i+1} = u;
    end;
    
     x = linspace(0,1,size(u,2));
     dt = 0.5/(64*64);
     Nf = size(u,1); % number of frames
     figure();
     axis([0 1 min([u(1,:) -1]) max([1 u(1,:)])])
     xlabel('x','Fontsize',14)
     ylabel('u(x)','Fontsize',14)
     hold on
     j = 1;
     title(sprintf('g_%d with %s', i+1, method_name));
     for i=1:Nf
         cla;
         for k = 1:3
            subplot(1, 3, k);
            u = ulist{k};
            
            plot(x,u(i,:));
            title(sprintf('Time = %0.4f', i*dt));
         end
         pause(0.000001)
         M(j)=getframe(gcf);
         j=j+1;
      end
      close('all');
    end
