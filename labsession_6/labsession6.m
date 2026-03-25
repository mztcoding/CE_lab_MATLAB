clc; clear; close all;

%% ── Physical Parameters ──────────────────────────────────
w        = 10;      % weight per unit length (N/m)
y0       = 5;       % minimum cable height (m)
x        = 50;      % horizontal distance (m)
y_target = 15;      % required height at x=50 (m)

%% ── Function Definitions ─────────────────────────────────
f  = @(TA) (TA./w).*cosh(w.*x./TA) + y0 - TA./w - y_target;
fp = @(TA) (1/w).*cosh(w.*x./TA) - (x./TA).*sinh(w.*x./TA) - 1/w;

%% ── Plot f(T_A) to identify root ─────────────────────────
TA_range = linspace(50, 2000, 2000);
figure;
plot(TA_range, f(TA_range), 'b-', 'LineWidth', 2); hold on;
yline(0, 'k--', 'LineWidth', 1);
xlabel('T_A (N/m)'); ylabel('f(T_A)');
title('Nonlinear Function f(T_A) vs Horizontal Tension');
legend('f(T_A)', 'f = 0'); grid on;

%% ── IQI Method ───────────────────────────────────────────
tol   = 0.0001;   % stopping tolerance (%)
maxit = 100;

[r_iqi, roots_iqi, ea_iqi, n_iqi] = iqi(f, 200, 350, 500, tol, maxit);
fprintf('IQI  :  T_A = %.6f N/m  (%d iterations)\n', r_iqi, n_iqi);

%% ── Newton–Raphson Method ────────────────────────────────
[r_nr, roots_nr, ea_nr, n_nr] = newton(f, fp, 350, tol, maxit);
fprintf('N-R  :  T_A = %.6f N/m  (%d iterations)\n', r_nr, n_nr);

%% ── Iteration Tables ─────────────────────────────────────
fprintf('\n--- IQI Iteration Table ---\n');
fprintf('%5s | %14s | %12s\n', 'Iter', 'T_A (N/m)', 'ea (%)');
for i = 1:length(roots_iqi)
    if isnan(ea_iqi(i))
        fprintf('%5d | %14.6f | %12s\n', i-1, roots_iqi(i), '---');
    else
        fprintf('%5d | %14.6f | %12.6f\n', i-1, roots_iqi(i), ea_iqi(i));
    end
end

fprintf('\n--- Newton-Raphson Iteration Table ---\n');
fprintf('%5s | %14s | %12s\n', 'Iter', 'T_A (N/m)', 'ea (%)');
for i = 1:length(roots_nr)
    if isnan(ea_nr(i))
        fprintf('%5d | %14.6f | %12s\n', i-1, roots_nr(i), '---');
    else
        fprintf('%5d | %14.6f | %12.6f\n', i-1, roots_nr(i), ea_nr(i));
    end
end

%% ── Catenary Cable Profile ───────────────────────────────
x_range = -50:0.5:100;
TA_sol  = r_iqi;
y_cable = (TA_sol/w).*cosh(w.*x_range./TA_sol) + y0 - TA_sol/w;

figure;
plot(x_range, y_cable, 'b-', 'LineWidth', 2.5); hold on;
plot(50, 15,  'rv', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
plot(0,  y0,  'gs', 'MarkerSize', 11, 'MarkerFaceColor', 'g');
yline(15, 'r--', 'LineWidth', 0.8);
xline(50, 'r--', 'LineWidth', 0.8);
xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Catenary Cable Profile  (T_A^* = %.2f N/m)', TA_sol));
legend('Cable profile', 'Target (50, 15)', 'Min point (0, 5)');
grid on;

%% ── Error Convergence Plot ───────────────────────────────
ea_iqi_plot = ea_iqi(~isnan(ea_iqi));
ea_nr_plot  = ea_nr(~isnan(ea_nr));

figure;
semilogy(1:length(ea_iqi_plot), ea_iqi_plot, 'bo-', 'LineWidth', 2, ...
         'MarkerSize', 8, 'DisplayName', 'IQI'); hold on;
semilogy(1:length(ea_nr_plot),  ea_nr_plot,  'rs--','LineWidth', 2, ...
         'MarkerSize', 8, 'DisplayName', 'Newton-Raphson');
yline(tol, 'k:', 'LineWidth', 1.2, 'DisplayName', '\epsilon_s = 0.0001%');
xlabel('Iteration'); ylabel('Approximate \epsilon_a (%)');
title('Error Convergence: IQI vs Newton-Raphson');
legend; grid on;

%% ── Iteration History Plot ───────────────────────────────
figure;
plot(0:length(roots_iqi)-1, roots_iqi, 'bo-', 'LineWidth', 2, ...
     'MarkerSize', 8, 'DisplayName', 'IQI'); hold on;
plot(0:length(roots_nr)-1,  roots_nr,  'rs--','LineWidth', 2, ...
     'MarkerSize', 8, 'DisplayName', 'Newton-Raphson');
yline(r_nr, 'k--', 'LineWidth', 1.2, 'DisplayName', 'True root T_A^*');
xlabel('Iteration index'); ylabel('T_A estimate (N/m)');
title('Iteration History: T_A Estimates Converging to Root');
legend; grid on;


%% =========================================================
%                     FUNCTIONS
%% =========================================================

function [root, roots_arr, ea_arr, n] = iqi(f, x0, x1, x2, tol, maxiter)
% IQI  Inverse Quadratic Interpolation.
%      Three initial guesses required (open method, no bracket needed).
%      Uses sliding-window update to advance the triple each iteration.
    roots_arr = [x0, x1, x2];
    ea_arr    = [NaN, NaN, NaN];
    n  = 0;
    xa = x0;  xb = x1;  xc = x2;

    for i = 1:maxiter
        n  = i;
        fa = f(xa);  fb = f(xb);  fc = f(xc);

        % Denominator guard – failure condition
        denom = (fb - fa) * (fc - fa) * (fc - fb);
        if abs(denom) < 1e-30
            error('IQI:DenomZero', ...
                  ['Denominator ≈ 0: function values nearly equal.\n' ...
                   'Restart with three better-separated initial guesses.']);
        end

        % IQI formula – Lagrange inverse parabola evaluated at y = 0
        xd = xa*(fb*fc) / ((fa-fb)*(fa-fc)) ...
           + xb*(fa*fc) / ((fb-fa)*(fb-fc)) ...
           + xc*(fa*fb) / ((fc-fa)*(fc-fb));

        ea = abs((xd - xc) / xd) * 100;
        roots_arr(end+1) = xd;
        ea_arr(end+1)    = ea;

        if f(xd) == 0 || ea < tol, break; end

        % Sliding-window update: drop oldest, add newest
        xa = xb;  xb = xc;  xc = xd;
    end
    root = xd;
end

% ---------------------------------------------------------

function [root, roots_arr, ea_arr, n] = newton(f, fp, x0, tol, maxiter)
% NEWTON  Newton-Raphson method (for convergence comparison).
    roots_arr = x0;
    ea_arr    = NaN;
    n  = 0;
    x  = x0;

    for i = 1:maxiter
        n  = i;
        xn = x - f(x)/fp(x);
        ea = abs((xn - x) / xn) * 100;
        roots_arr(end+1) = xn;
        ea_arr(end+1)    = ea;
        if ea < tol, break; end
        x = xn;
    end
    root = xn;
end
