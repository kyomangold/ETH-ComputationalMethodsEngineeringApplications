errors = load("errors_fd.txt");
resolutions = load("resolutions_fd.txt");
rate = polyfit(log(resolutions), log(errors), 1);
rate = rate(1);
disp("Approximate order: ");
disp(round(rate, 2));
N = length(errors);

loglog(resolutions, errors, '-o')
hold on
loglog(resolutions, errors(N-1)*resolutions.^rate/(resolutions(N-1).^rate), '--')
xlabel('# cells')
ylabel('L^1 error')
axis([resolutions(1) resolutions(N) -inf inf]);
label = strcat('N^{', num2str(rate,3), '}');
legend('||u_N-u||_\infty', label, 'Location', 'northeast')
grid on
