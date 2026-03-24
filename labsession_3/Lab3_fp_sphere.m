% =========================================================
%  LAB 3 – Problem 3: False Position – Floating Sphere
%  Finds the above-water height h of a sphere floating
%  in equilibrium using Archimedes' principle.
%  r = 1 m, rho_s = 200 kg/m^3, rho_w = 1000 kg/m^3
%  f(h) = (pi*h^2/3)*(3r-h) - (4*pi*r^3/3)*(rho_s/rho_w)
%
%  REQUIRES: false_position.m in the same folder.
% =========================================================
clc; clear; close all;

r     = 1;      % sphere radius (m)
rho_s = 200;    % sphere density (kg/m^3)
rho_w = 1000;   % water density  (kg/m^3)

f = @(h) (pi.*h.^2/3).*(3*r - h) - (4*pi*r^3/3)*(rho_s/rho_w);

xl = 0;  xu = 2*r;   % h must be in (0, 2r)
tol = 0.001;  max_iter = 100;

[xr, roots_arr, ea_arr, n_iter] = false_position(f, xl, xu, tol, max_iter);
fprintf('Above-water height h = %.6f m  (%d iterations)\n', xr, n_iter);
fprintf('Submerged depth      = %.6f m\n', 2*r - xr);

%% ===== PLOT 1: Function graph with markers =====
h_plot = linspace(0, 2, 400);
figure;
plot(h_plot, f(h_plot), 'b-', 'LineWidth',2);  hold on;
yline(0,'k--');
for k = 1:length(roots_arr)
    plot(roots_arr(k), f(roots_arr(k)), 'ro', ...
         'MarkerFaceColor','r','MarkerSize',7);
    text(roots_arr(k), f(roots_arr(k))+0.01, num2str(k), ...
         'FontSize',8,'Color','r','HorizontalAlignment','center');
end
xlabel('h (m)');  ylabel('f(h)');
title('False Position – Floating Sphere (Problem 5.22)');  grid on;

%% ===== PLOT 2: Convergence =====
figure;
semilogy(2:length(ea_arr), ea_arr(2:end), 'd-', 'LineWidth',1.5, ...
         'MarkerFaceColor','b');
xlabel('Iteration');  ylabel('Approx. % Relative Error');
title('False Position Convergence – Floating Sphere');  grid on;
