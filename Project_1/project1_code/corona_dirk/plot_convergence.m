clc
close all
clear all

FolderName = 'build/';

error = load([FolderName,'errors.txt']);
numbers = load([FolderName,'numbers.txt']);
times = load([FolderName,'walltimes.txt']);

if (max(error) == 0 || max(numbers) == 0)
    disp("errors.txt or numbers.txt is all zeros!");
    return
end

N = length(error);
rate = polyfit(log(numbers(3:N)), log(error(3:N)), 1);
rate = rate(1);
disp("Error scales with number of cells with order: ");
disp(rate);
loglog(numbers, error, '-o')
hold on
loglog(numbers, error(N-1)*numbers.^rate/(numbers(N-1).^rate), '--')
xlabel('# cells')
ylabel('L^1 error')
axis([numbers(1) numbers(N) -inf inf]);
label = strcat('N^{', num2str(rate,3), '}');
legend('|u(T)-U_N|', label, 'Location', 'northeast','FontSize',13)
grid on

if (max(times) > 0)
    disp('Press any key to continue...')
    pause
    cla

    rate = polyfit(log(numbers(3:N)), log(times(3:N)), 1);
    rate = rate(1);
    disp("Runtime scales with number of cells with order: ");
    disp(rate);
    loglog(numbers, times, '-o')
    hold on
    loglog(numbers, times(N-1)*numbers.^rate/(numbers(N-1).^rate), '--')
    xlabel('# cells')
    ylabel('ns')
    axis([numbers(1) numbers(N) -inf inf]);
    label = strcat('N^{', num2str(rate,3), '}');
    legend('Walltime', label, 'Location', 'northwest','FontSize',13)
    grid on

    disp('Press any key to continue...')
    pause
    cla

    rate = polyfit(log(times(3:N)), log(error(3:N)), 1);
    rate = rate(1);
    disp("Error scales with runtime with order: ");
    disp(rate);
    loglog(times, error, '-o')
    hold on
    loglog(times, error(N-1)*times.^rate/(times(N-1).^rate), '--')
    xlabel('ns')
    ylabel('L^1 error')
    axis([times(1) times(N) -inf inf]);
    label = strcat('t^{', num2str(rate,3), '}');
    legend('error', label, 'Location', 'northeast','FontSize',13)
    grid on
end

