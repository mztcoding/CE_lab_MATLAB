% =========================================================
%  LAB 3 – Problem 1: False Position vs Bisection
%  Finds the root of f(x) = x^10 - 1 on [0, 1.3].
%  Demonstrates stagnation of the false position method
%  on highly curved functions.
%
%  REQUIRES: false_position.m in the same folder.
% =========================================================
clc; clear; close all;

f    = @(x) x.^10 - 1;
xl   = 0;   xu = 1.3;
tol  = 0.001;   max_iter = 100;

%% ===== FALSE POSITION =====
[xr_fp, roots_fp, ea_fp, n_fp] = false_position(f, xl, xu, tol, max_iter);
fprintf('False Position: root = %.8f  (%d iterations)\n', xr_fp, n_fp);

%% ===== BISECTION =====
xr_b = xl;  xl_b = xl;  xu_b = xu;
roots_bs = [];  ea_bs = [];
for i = 1:max_iter
    xr_new = (xl_b + xu_b) / 2;
    if i > 1
        ea_b = abs((xr_new - xr_b) / xr_new) * 100;
    else
        ea_b = inf;
    end
    roots_bs(end+1) = xr_new;
    ea_bs(end+1)    = ea_b;
    if f(xl_b)*f(xr_new) < 0,  xu_b = xr_new;
    else,                       xl_b = xr_new;
    end
    xr_b = xr_new;
    if ea_b < tol,  break;  end
end
fprintf('Bisection:       root = %.8f  (%d iterations)\n', xr_b, i);

%% ===== PLOT 1: Function with FP markers =====
x_plot = linspace(0, 1.3, 500);
figure;
plot(x_plot, f(x_plot), 'b-', 'LineWidth', 2);  hold on;
yline(0,'k--');
for k = 1:length(roots_fp)
    plot(roots_fp(k), f(roots_fp(k)), 'ro', 'MarkerFaceColor','r','MarkerSize',6);
    text(roots_fp(k), f(roots_fp(k))+0.03, num2str(k), 'FontSize',7,'Color','r', ...
         'HorizontalAlignment','center');
end
xlabel('x');  ylabel('f(x) = x^{10} - 1');
title('False Position – f(x) = x^{10} - 1');  grid on;

%% ===== PLOT 2: Comparative convergence (error vs iteration) =====
figure;
semilogy(2:length(ea_fp), ea_fp(2:end), 'r-o', 'LineWidth',1.5, ...
         'DisplayName','False Position');  hold on;
semilogy(2:length(ea_bs), ea_bs(2:end), 'b-s', 'LineWidth',1.5, ...
         'DisplayName','Bisection');
xlabel('Iteration');  ylabel('Approx. % Relative Error');
title('Convergence Comparison: False Position vs Bisection');
legend;  grid on;

%% ===== PLOT 3: Successive approximations on same graph =====
figure;
plot(x_plot, f(x_plot), 'k-', 'LineWidth', 2);  hold on;
yline(0,'k--');
plot(roots_fp, f(roots_fp), 'r^', 'MarkerFaceColor','r','MarkerSize',7, ...
     'DisplayName','False Position');
plot(roots_bs, f(roots_bs), 'bs', 'MarkerFaceColor','b','MarkerSize',7, ...
     'DisplayName','Bisection');
xlabel('x');  ylabel('f(x)');
title('Root Approximations: FP vs Bisection on x^{10}-1');
legend;  grid on;
