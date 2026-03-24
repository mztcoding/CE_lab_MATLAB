% =========================================================
%  LAB 3 – Problem 2: False Position – Rocket Velocity
%  Finds the time t at which a rocket reaches v = 750 m/s.
%  v(t) = u*ln(m0/(m0-q*t)) - g*t
%  u=1800 m/s, m0=160000 kg, q=2600 kg/s, g=9.81 m/s^2
%  Bracket: [10, 50] s   Tolerance: 1%
%
%  REQUIRES: false_position.m in the same folder.
% =========================================================
clc; clear; close all;

u  = 1800;  m0 = 160000;  q = 2600;  g = 9.81;  v_target = 750;
f  = @(t) u.*log(m0./(m0 - q.*t)) - g.*t - v_target;

xl = 10;   xu = 50;
tol = 1;   max_iter = 100;   % within 1% of true value

[xr, roots_arr, ea_arr, n_iter] = false_position(f, xl, xu, tol, max_iter);
fprintf('Root: t = %.6f s  (after %d iterations)\n', xr, n_iter);
fprintf('Verification: v(%.4f) = %.4f m/s\n', xr, ...
        u*log(m0/(m0-q*xr)) - g*xr);

%% ===== PLOT 1: Function graph with markers =====
t_plot = linspace(10, 50, 500);
figure;
plot(t_plot, f(t_plot), 'b-', 'LineWidth',2);  hold on;
yline(0,'k--');
for k = 1:length(roots_arr)
    plot(roots_arr(k), f(roots_arr(k)), 'ro', ...
         'MarkerFaceColor','r','MarkerSize',7);
    text(roots_arr(k), f(roots_arr(k))+3, num2str(k), ...
         'FontSize',8,'Color','r','HorizontalAlignment','center');
end
xlabel('t (s)');  ylabel('f(t) = v(t) - 750');
title('False Position – Rocket Velocity Problem');  grid on;

%% ===== PLOT 2: Convergence =====
figure;
semilogy(2:length(ea_arr), ea_arr(2:end), 's-', 'LineWidth',1.5, ...
         'MarkerFaceColor','b');
xlabel('Iteration');  ylabel('Approx. % Relative Error');
title('False Position Convergence – Rocket Problem');  grid on;
