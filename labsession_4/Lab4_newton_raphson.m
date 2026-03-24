% =========================================================
%  LAB 4 – Newton-Raphson Method
%  Spring-mass deflection problem.
%  k1=40000 N/m, k2=40 N/m^1.5, m=0.095 kg, g=9.81, h=0.43 m
%  f(d)  = (2k2/5)*d^2.5 + 0.5*k1*d^2 - mg*d - mgh
%  f'(d) = k2*d^1.5 + k1*d - mg
%
%  Enter an initial guess when prompted (try 0.01).
% =========================================================
clc;
clear;

% ===== GIVEN PARAMETERS =====
k1  = 40000;
k2  = 40;
m   = 0.095;
g   = 9.81;
h   = 0.43;
mg  = m * g;
mgh = mg * h;

% Function and its derivative
f  = @(d) (2*k2/5)*d.^(2.5) + 0.5*k1*d.^2 - mg*d - mgh;
df = @(d) k2*d.^(1.5) + k1*d - mg;

% ===== USER INPUT =====
d        = input('Enter initial guess: ');
tol      = 1e-6;
max_iter = 100;
errors   = zeros(1, max_iter);

fprintf('\nNewton-Raphson Method\n');

% ===== ITERATION =====
for i = 1:max_iter
    d_new      = d - f(d)/df(d);
    errors(i)  = abs(d_new - d);

    fprintf('Iteration %d: %.8f\n', i, d_new);

    if errors(i) < tol
        break;
    end
    d = d_new;
end

root = d_new;
fprintf('\nFinal Root = %.8f m\n', root);

% ===== PLOT 1: Convergence =====
figure;
plot(1:i, errors(1:i), '-o');
xlabel('Iteration');  ylabel('Error');
title('Newton-Raphson Convergence');
grid on;

% ===== PLOT 2: Function with tangent line at root =====
d_vals  = linspace(0.001, 0.008, 500);
tangent = f(root) + df(root).*(d_vals - root);
figure;
plot(d_vals, f(d_vals), 'LineWidth', 1.5);  hold on;
plot(d_vals, tangent, 'r--', 'LineWidth', 1.5);
yline(0, 'k--');
xlabel('d (m)');  ylabel('f(d)');
legend('f(d)', 'Tangent Line', 'Location', 'Best');
title('Newton-Raphson Tangent at Root');
grid on;
