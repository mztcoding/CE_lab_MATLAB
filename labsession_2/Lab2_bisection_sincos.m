% =========================================================
%  LAB 2 – Problem 2: Bisection Method
%  Finds a root of f(x) = sin(10x) + cos(3x).
%  The function is plotted first so the user can identify
%  a valid sign-change bracket before entering xl and xu.
%  Stopping criterion: ea < 0.01%
% =========================================================
clc; clear; close all;

%% ===== PROBLEM DEFINITION =====
f   = @(x) sin(10*x) + cos(3*x);
tol = 0.0001;
max_iter = 100;

%% ===== PLOT FUNCTION TO IDENTIFY BRACKET =====
x_plot = linspace(0, 5, 1000);
figure;
plot(x_plot, f(x_plot), 'b-', 'LineWidth', 2);  hold on;
yline(0, 'k--');
xlabel('x');  ylabel('f(x)');
title('f(x) = sin(10x) + cos(3x)  –  Identify a sign-change interval');
grid on;

%% ===== USER INPUT FOR BRACKET =====
xl = input('Enter lower bracket xl: ');
xu = input('Enter upper bracket xu: ');

if f(xl)*f(xu) >= 0
    error('No sign change in the given interval.');
end

%% ===== BISECTION ITERATION =====
xr_old    = xl;
ea        = inf;
roots_arr = [];
ea_arr    = [];
iter      = 0;

while ea > tol && iter < max_iter
    iter = iter + 1;
    xr   = (xl + xu) / 2;

    if iter > 1
        ea = abs((xr - xr_old) / xr) * 100;
    end

    roots_arr(end+1) = xr;
    ea_arr(end+1)    = ea;

    fprintf('Iter %2d: xr = %.6f,  ea = %.4f%%\n', iter, xr, ea);

    if f(xl)*f(xr) < 0
        xu = xr;
    elseif f(xr)*f(xu) < 0
        xl = xr;
    else
        ea = 0;  break;
    end
    xr_old = xr;
end

fprintf('\nRoot: x = %.6f  (after %d iterations)\n', xr, iter);

%% ===== OVERLAY ITERATION POINTS ON EXISTING FIGURE =====
figure(1);
for k = 1:length(roots_arr)
    plot(roots_arr(k), f(roots_arr(k)), 'ro', ...
         'MarkerFaceColor', 'r', 'MarkerSize', 6);
    text(roots_arr(k), f(roots_arr(k)) + 0.06, num2str(k), ...
         'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', 'r');
end
title('Bisection on f(x) = sin(10x) + cos(3x) – Iteration markers');

%% ===== CONVERGENCE PLOT =====
figure;
semilogy(2:length(ea_arr), ea_arr(2:end), 's-', 'LineWidth', 1.5, ...
         'MarkerFaceColor', 'b');
xlabel('Iteration');  ylabel('Approximate % Relative Error');
title('Bisection Convergence – sin(10x)+cos(3x)');
grid on;
