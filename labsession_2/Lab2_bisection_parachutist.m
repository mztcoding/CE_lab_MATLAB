% =========================================================
%  LAB 2 – Problem 1: Bisection Method
%  Finds the drag coefficient c for a parachutist of
%  mass m = 68.1 kg to reach v = 40 m/s after t = 10 s.
%  User enters the bracketing interval [xl, xu].
%  Stopping criterion: ea < 0.01%
% =========================================================
clc; clear; close all;

%% ===== PROBLEM DEFINITION =====
m = 68.1;  g = 9.8;  v_target = 40;  t = 10;
f = @(c) (g*m./c) .* (1 - exp(-(c/m)*t)) - v_target;

tol      = 0.0001;   % 0.01% stopping criterion
max_iter = 100;

%% ===== USER INPUT FOR BRACKET =====
xl = input('Enter lower bracket xl: ');
xu = input('Enter upper bracket xu: ');

if f(xl)*f(xu) >= 0
    error('f(xl) and f(xu) must have opposite signs.');
end

%% ===== BISECTION ITERATION =====
xr_old    = xl;
ea        = inf;
roots_arr = [];
ea_arr    = [];
iter      = 0;

while ea > tol && iter < max_iter
    iter  = iter + 1;
    xr    = (xl + xu) / 2;

    if iter > 1
        ea = abs((xr - xr_old) / xr) * 100;
    end

    roots_arr(end+1) = xr;
    ea_arr(end+1)    = ea;

    fprintf('Iter %2d: xr = %.6f,  ea = %.4f%%\n', iter, xr, ea);

    test = f(xl) * f(xr);
    if test < 0
        xu = xr;
    elseif test > 0
        xl = xr;
    else
        ea = 0;  break;
    end
    xr_old = xr;
end

fprintf('\nRoot found: c = %.6f  (after %d iterations)\n', xr, iter);

%% ===== PLOT 1: Function with iteration markers =====
c_vals = linspace(4, 20, 400);
figure;
plot(c_vals, f(c_vals), 'b-', 'LineWidth', 2);  hold on;
yline(0, 'k--', 'LineWidth', 1);
for k = 1:length(roots_arr)
    plot(roots_arr(k), f(roots_arr(k)), 'ro', 'MarkerFaceColor', 'r', ...
         'MarkerSize', 7);
    text(roots_arr(k), f(roots_arr(k)) + 0.8, num2str(k), ...
         'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', 'r');
end
xlabel('c  (kg/s)');  ylabel('f(c)');
title('Bisection – Parachutist Drag Coefficient');
grid on;

%% ===== PLOT 2: Approximate % relative error vs iteration =====
figure;
semilogy(2:length(ea_arr), ea_arr(2:end), 's-', 'LineWidth', 1.5, ...
         'MarkerFaceColor', 'b');
xlabel('Iteration');  ylabel('Approximate % Relative Error');
title('Bisection Convergence – Parachutist Problem');
grid on;
