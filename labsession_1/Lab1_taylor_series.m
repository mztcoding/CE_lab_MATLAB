% =========================================================
%  LAB 1 – Problem 3: Taylor Series Approximation of sin(x)
%  Plots the Maclaurin series for N = 1 to 10 terms over
%  [-360, 360] degrees, analyses convergence at x = 60 deg,
%  and displays an error summary table.
% =========================================================
clc; clear; close all;

%% ===== PARAMETERS =====
Nterms  = 10;
x_deg   = -360:1:360;
x       = deg2rad(x_deg);

%% ===== TRUE FUNCTION =====
f_true  = sin(x);

%% ===== MACLAURIN SERIES APPROXIMATIONS =====
T = zeros(Nterms, length(x));
for N = 1:Nterms
    Tsum = 0;
    for n = 0:N-1
        Tsum = Tsum + ((-1)^n / factorial(2*n+1)) .* x.^(2*n+1);
    end
    T(N,:) = Tsum;
end

%% ===== PLOT 1: sin(x) and all Taylor approximations =====
figure;
plot(x_deg, f_true, 'k', 'LineWidth', 3, 'DisplayName', 'True sin(x)');
hold on;
for N = 1:Nterms
    plot(x_deg, T(N,:), '--', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%d terms', N));
end
xlabel('x (degrees)');  ylabel('sin(x)');
title('sin(x) and Maclaurin Series Approximations');
legend('Location','best');
grid on;  xlim([-360 360]);  ylim([-1.5 1.5]);

%% ===== VALUE AT x = 60 degrees =====
x60       = deg2rad(60);
true_val  = sin(x60);
approx_val = zeros(1, Nterms);
for N = 1:Nterms
    for n = 0:N-1
        approx_val(N) = approx_val(N) + ...
            ((-1)^n / factorial(2*n+1)) * x60^(2*n+1);
    end
end

%% ===== PLOT 2: Convergence at x = 60 degrees =====
figure;
plot(1:Nterms, approx_val, 'o-', 'LineWidth', 1.5);  hold on;
yline(true_val, 'r--', 'LineWidth', 2);
xlabel('Number of Terms');  ylabel('sin(60°)');
title('Convergence of Taylor Series at x = 60°');
legend('Taylor Approximation','True Value');
grid on;

%% ===== ERROR ANALYSIS =====
true_error = true_val - approx_val;
abs_error  = abs(true_error);
rel_error  = abs_error ./ abs(true_val);

%% ===== PLOT 3: Absolute Error =====
figure;
semilogy(1:Nterms, abs_error, 's-', 'LineWidth', 1.5);
xlabel('Number of Terms');  ylabel('Absolute Error');
title('Absolute Error vs. Number of Terms');
grid on;

%% ===== PLOT 4: Relative Error =====
figure;
semilogy(1:Nterms, rel_error, 'd-', 'LineWidth', 1.5);
xlabel('Number of Terms');  ylabel('Relative Error');
title('Relative Error vs. Number of Terms');
grid on;

%% ===== DISPLAY TABLE =====
fprintf('\nError Summary at x = 60 degrees\n');
fprintf('%-8s  %-14s  %-14s  %-14s\n', ...
        'Terms','Approx Value','Abs Error','Rel Error');
for N = 1:Nterms
    fprintf('%-8d  %-14.8f  %-14.2e  %-14.2e\n', ...
            N, approx_val(N), abs_error(N), rel_error(N));
end
