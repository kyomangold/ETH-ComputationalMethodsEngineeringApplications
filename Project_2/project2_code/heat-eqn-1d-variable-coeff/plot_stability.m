function plot_stability(method)

    coefficients = {'a_1', 'a_2', 'a_3'};

    N = 127;
    dxSquared = 1.0 / (N+1)^2;

    timesteps = 0.5*dxSquared * [128., 8., 1.];
    

    for i = 1:3
        
        for j = 1:3
            filename = sprintf('%s_stability_u_a%d_dt%d.txt',method, i-1, j-1);

            u = load(filename);
            x = linspace(0, 1, size(u, 2));

            subplot(3, 3, 3*(i-1) + j);
            plot(x, u(end,:));
            title(sprintf('a_{%d}, \\Delta t_{%d}', i, j));
        end;
    end;
    
       
