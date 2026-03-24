%% ============================================================
%  LAB SESSION 5 — Secant Method and Modified Secant Method
%  Projectile Trajectory Problem (Chapra Problem 6.21)
%
%  Finds launch angles for a ball thrown to a catcher at:
%    x = 90 m,  y_c = 1.0 m,  v0 = 30 m/s,  y0 = 1.8 m
%
%  Methods implemented:
%    1. Secant Method
%    2. Modified Secant Method
%    3. Newton-Raphson (for comparison)
%
%  All figures saved to graphs/ folder as PDF.
%% ============================================================
clc; clear; close all;

% Create output folder if it doesn't exist
if ~exist('graphs', 'dir')
    mkdir('graphs');
end

%% ── 1. PHYSICAL PARAMETERS ──────────────────────────────────────────────
g     = 9.81;      % gravitational acceleration (m/s^2)
v0    = 30;        % initial speed (m/s)
xc    = 90;        % horizontal distance to catcher (m)  [renamed from x to avoid clash]
y0    = 1.8;       % release height (m)
yc    = 1.0;       % catcher height (m)
tol   = 0.0001;    % stopping tolerance (%)
maxit = 100;       % maximum iterations

%% ── 2. FUNCTION DEFINITIONS ─────────────────────────────────────────────
% Nonlinear equation  f(theta) = 0
f  = @(th) tan(th).*xc ...
         - (g ./ (2*v0^2 .* cos(th).^2)) .* xc^2 ...
         + y0 - yc;

% Analytical derivative  f'(theta)  — used only by Newton-Raphson
fp = @(th) xc ./ cos(th).^2 ...
         - (g .* xc^2 .* sin(th)) ./ (v0^2 .* cos(th).^3);

%% ── 3. FUNCTION PLOT ────────────────────────────────────────────────────
th_deg = linspace(1, 89, 1000);
th_rad = deg2rad(th_deg);
fvals  = f(th_rad);

figure('Name','f(theta) plot','NumberTitle','off');
plot(th_deg, fvals, 'b-', 'LineWidth', 2); hold on;
yline(0, 'k--', 'LineWidth', 0.9);
xline(38, 'r:', 'LineWidth', 1.2);
xline(52, 'g:', 'LineWidth', 1.2);
xlabel('\theta_0 (degrees)');
ylabel('f(\theta_0)');
title('Nonlinear Function f(\theta_0) vs Launch Angle');
legend('f(\theta_0)', 'f = 0', 'Root 1 \approx 38°', 'Root 2 \approx 52°', ...
       'Location','northwest');
grid on;
xlim([0 90]); ylim([-60 60]);
saveas(gcf, 'graphs/lab5_fig1_fplot.pdf');

%% ── 4. SECANT METHOD ────────────────────────────────────────────────────
fprintf('\n===== SECANT METHOD =====\n');

[r1_s, h1_s, e1_s, n1_s] = secant_method(f, 20, 30, tol, maxit);
fprintf('Root 1: theta = %.6f deg   (%d iterations)\n', r1_s, n1_s);

[r2_s, h2_s, e2_s, n2_s] = secant_method(f, 60, 70, tol, maxit);
fprintf('Root 2: theta = %.6f deg   (%d iterations)\n', r2_s, n2_s);

fprintf('\n--- Secant Method: Root 1 ---\n');
fprintf('%6s  %14s  %14s\n','Iter','theta (deg)','ea (%)');
for k = 1:length(h1_s)
    if isnan(e1_s(k))
        fprintf('%6d  %14.6f  %14s\n', k-1, h1_s(k), '---');
    else
        fprintf('%6d  %14.6f  %14.6f\n', k-1, h1_s(k), e1_s(k));
    end
end

fprintf('\n--- Secant Method: Root 2 ---\n');
fprintf('%6s  %14s  %14s\n','Iter','theta (deg)','ea (%)');
for k = 1:length(h2_s)
    if isnan(e2_s(k))
        fprintf('%6d  %14.6f  %14s\n', k-1, h2_s(k), '---');
    else
        fprintf('%6d  %14.6f  %14.6f\n', k-1, h2_s(k), e2_s(k));
    end
end

%% ── 5. MODIFIED SECANT METHOD ───────────────────────────────────────────
fprintf('\n===== MODIFIED SECANT METHOD =====\n');
delta = 1e-4;

[r1_ms, h1_ms, e1_ms, n1_ms] = mod_secant_method(f, 30, delta, tol, maxit);
fprintf('Root 1: theta = %.6f deg   (%d iterations)\n', r1_ms, n1_ms);

[r2_ms, h2_ms, e2_ms, n2_ms] = mod_secant_method(f, 65, delta, tol, maxit);
fprintf('Root 2: theta = %.6f deg   (%d iterations)\n', r2_ms, n2_ms);

fprintf('\n--- Modified Secant Method: Root 1 ---\n');
fprintf('%6s  %14s  %14s\n','Iter','theta (deg)','ea (%)');
for k = 1:length(h1_ms)
    if isnan(e1_ms(k))
        fprintf('%6d  %14.6f  %14s\n', k-1, h1_ms(k), '---');
    else
        fprintf('%6d  %14.6f  %14.6f\n', k-1, h1_ms(k), e1_ms(k));
    end
end

fprintf('\n--- Modified Secant Method: Root 2 ---\n');
fprintf('%6s  %14s  %14s\n','Iter','theta (deg)','ea (%)');
for k = 1:length(h2_ms)
    if isnan(e2_ms(k))
        fprintf('%6d  %14.6f  %14s\n', k-1, h2_ms(k), '---');
    else
        fprintf('%6d  %14.6f  %14.6f\n', k-1, h2_ms(k), e2_ms(k));
    end
end

%% ── 6. NEWTON-RAPHSON METHOD ────────────────────────────────────────────
fprintf('\n===== NEWTON-RAPHSON METHOD =====\n');

[r1_nr, h1_nr, e1_nr, n1_nr] = newton_raphson(f, fp, 30, tol, maxit);
fprintf('Root 1: theta = %.6f deg   (%d iterations)\n', r1_nr, n1_nr);

[r2_nr, h2_nr, e2_nr, n2_nr] = newton_raphson(f, fp, 65, tol, maxit);
fprintf('Root 2: theta = %.6f deg   (%d iterations)\n', r2_nr, n2_nr);

%% ── 7. ERROR CONVERGENCE PLOT — ROOT 1 ──────────────────────────────────
ea_s1  = e1_s(~isnan(e1_s));
ea_ms1 = e1_ms(~isnan(e1_ms));
ea_nr1 = e1_nr(~isnan(e1_nr));

figure('Name','Error Convergence Root 1','NumberTitle','off');
semilogy(1:length(ea_s1),  ea_s1,  'bo-',  'LineWidth',1.8, 'MarkerSize',6); hold on;
semilogy(1:length(ea_ms1), ea_ms1, 'rs--', 'LineWidth',1.8, 'MarkerSize',6);
semilogy(1:length(ea_nr1), ea_nr1, 'g^:',  'LineWidth',1.8, 'MarkerSize',6);
yline(tol, 'k-.', 'LineWidth', 1.0);                      % tolerance line
xlabel('Iteration');
ylabel('Approximate \varepsilon_a (%)');
title('Error Convergence — Root 1 (\theta^* \approx 37.96°)');
legend('Secant','Modified Secant','Newton-Raphson', ...
       ['\varepsilon_s = ' num2str(tol) '%'], 'Location','northeast');
grid on;
saveas(gcf, 'graphs/lab5_fig2_error_r1.pdf');

%% ── 8. ERROR CONVERGENCE PLOT — ROOT 2 ──────────────────────────────────
ea_s2  = e2_s(~isnan(e2_s));
ea_ms2 = e2_ms(~isnan(e2_ms));
ea_nr2 = e2_nr(~isnan(e2_nr));

figure('Name','Error Convergence Root 2','NumberTitle','off');
semilogy(1:length(ea_s2),  ea_s2,  'bo-',  'LineWidth',1.8, 'MarkerSize',6); hold on;
semilogy(1:length(ea_ms2), ea_ms2, 'rs--', 'LineWidth',1.8, 'MarkerSize',6);
semilogy(1:length(ea_nr2), ea_nr2, 'g^:',  'LineWidth',1.8, 'MarkerSize',6);
yline(tol, 'k-.', 'LineWidth', 1.0);
xlabel('Iteration');
ylabel('Approximate \varepsilon_a (%)');
title('Error Convergence — Root 2 (\theta^* \approx 51.53°)');
legend('Secant','Modified Secant','Newton-Raphson', ...
       ['\varepsilon_s = ' num2str(tol) '%'], 'Location','northeast');
grid on;
saveas(gcf, 'graphs/lab5_fig3_error_r2.pdf');

%% ── 9. BALL TRAJECTORIES ────────────────────────────────────────────────
x_arr   = linspace(0, xc, 500);
y1_traj = traj(r1_s, x_arr, g, v0, y0);
y2_traj = traj(r2_s, x_arr, g, v0, y0);

figure('Name','Ball Trajectories','NumberTitle','off');
plot(x_arr, y1_traj, 'b-',  'LineWidth', 2.5); hold on;
plot(x_arr, y2_traj, 'r--', 'LineWidth', 2.5);
plot(0,   y0, 'ko', 'MarkerSize', 8,  'MarkerFaceColor','k');
plot(xc,  yc, 'g^', 'MarkerSize', 10, 'MarkerFaceColor','g');
xlabel('Horizontal distance x (m)');
ylabel('Height y (m)');
title('Ball Trajectories for Both Valid Launch Angles');
legend(sprintf('Root 1: \\theta_0 = %.2f° (low angle)',  r1_s), ...
       sprintf('Root 2: \\theta_0 = %.2f° (high angle)', r2_s), ...
       'Release point (0, 1.8 m)', ...
       'Catcher (90 m, 1.0 m)', ...
       'Location','northwest');
grid on;
xlim([-2 95]); ylim([0 35]);
saveas(gcf, 'graphs/lab5_fig4_trajectories.pdf');

%% ── 10. ITERATION HISTORY PLOT ───────────────────────────────────────────
figure('Name','Iteration History','NumberTitle','off');

subplot(2,1,1);
plot(0:length(h1_s)-1,  h1_s,  'bo-',  'LineWidth',1.6, 'MarkerSize',5); hold on;
plot(0:length(h1_ms)-1, h1_ms, 'rs--', 'LineWidth',1.6, 'MarkerSize',5);
plot(0:length(h1_nr)-1, h1_nr, 'g^:',  'LineWidth',1.6, 'MarkerSize',5);
yline(r1_s, 'k--', 'LineWidth', 1.0);
xlabel('Iteration'); ylabel('\theta_0 estimate (deg)');
title('Iteration History — Root 1 (\theta^* \approx 37.96°)');
legend('Secant','Modified Secant','Newton-Raphson','True root','Location','southeast');
grid on;

subplot(2,1,2);
plot(0:length(h2_s)-1,  h2_s,  'bo-',  'LineWidth',1.6, 'MarkerSize',5); hold on;
plot(0:length(h2_ms)-1, h2_ms, 'rs--', 'LineWidth',1.6, 'MarkerSize',5);
plot(0:length(h2_nr)-1, h2_nr, 'g^:',  'LineWidth',1.6, 'MarkerSize',5);
yline(r2_s, 'k--', 'LineWidth', 1.0);
xlabel('Iteration'); ylabel('\theta_0 estimate (deg)');
title('Iteration History — Root 2 (\theta^* \approx 51.53°)');
legend('Secant','Modified Secant','Newton-Raphson','True root','Location','southeast');
grid on;

saveas(gcf, 'graphs/lab5_fig5_iter_history.pdf');

%% ── 11. MODIFIED SECANT vs NEWTON-RAPHSON COMPARISON ────────────────────
figure('Name','ModSec vs NR','NumberTitle','off');

subplot(1,2,1);
semilogy(1:length(ea_ms1), ea_ms1, 'rs--', 'LineWidth',1.8, 'MarkerSize',6); hold on;
semilogy(1:length(ea_nr1), ea_nr1, 'g^:',  'LineWidth',1.8, 'MarkerSize',6);
xlabel('Iteration'); ylabel('\varepsilon_a (%)');
title('Root 1 (\theta^* \approx 37.96°)');
legend('Modified Secant','Newton-Raphson','Location','northeast');
grid on;

subplot(1,2,2);
semilogy(1:length(ea_ms2), ea_ms2, 'rs--', 'LineWidth',1.8, 'MarkerSize',6); hold on;
semilogy(1:length(ea_nr2), ea_nr2, 'g^:',  'LineWidth',1.8, 'MarkerSize',6);
xlabel('Iteration'); ylabel('\varepsilon_a (%)');
title('Root 2 (\theta^* \approx 51.53°)');
legend('Modified Secant','Newton-Raphson','Location','northeast');
grid on;

sgtitle('Modified Secant vs Newton-Raphson: Error Convergence');
saveas(gcf, 'graphs/lab5_fig6_modsec_vs_nr.pdf');

%% ── 12. SUMMARY TABLE ────────────────────────────────────────────────────
fprintf('\n===== COMPARATIVE SUMMARY =====\n');
fprintf('%-20s  %-10s  %-14s  %-14s  %-14s\n', ...
        'Method','Derivative','Inputs','Iters Root 1','Iters Root 2');
fprintf('%s\n', repmat('-',1,78));
fprintf('%-20s  %-10s  %-14s  %-14d  %-14d\n', ...
        'Secant',         'No',  'x_{i-1}, x_i', n1_s,  n2_s);
fprintf('%-20s  %-10s  %-14s  %-14d  %-14d\n', ...
        'Modified Secant','No',  'x_0, delta',   n1_ms, n2_ms);
fprintf('%-20s  %-10s  %-14s  %-14d  %-14d\n', ...
        'Newton-Raphson', 'Yes', 'x_0',          n1_nr, n2_nr);
fprintf('%s\n', repmat('-',1,78));
fprintf('\nAll figures saved to graphs/ folder.\n');

%% ════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════

% ── Secant Method ────────────────────────────────────────────────────────
function [root, hist, ea_hist, n] = secant_method(f, x0d, x1d, tol, maxiter)
% SECANT_METHOD  Secant method with two initial guesses (in degrees).
    x0 = deg2rad(x0d);
    x1 = deg2rad(x1d);
    hist    = [x0d, x1d];
    ea_hist = [NaN, NaN];
    n       = 0;
    x2      = x1;    % pre-assign so root is defined if maxiter=0

    if abs(f(x0) - f(x1)) < eps
        error('secant_method: f(x0) ≈ f(x1). Choose different initial guesses.');
    end

    for i = 1:maxiter
        n     = i;
        denom = f(x0) - f(x1);
        if abs(denom) < eps
            warning('secant_method: near-zero denominator at iteration %d. Stopping.', i);
            break;
        end
        x2 = x1 - f(x1) * (x0 - x1) / denom;
        ea = abs((x2 - x1) / x2) * 100;
        hist(end+1)    = rad2deg(x2);   %#ok<AGROW>
        ea_hist(end+1) = ea;            %#ok<AGROW>
        if f(x2) == 0 || ea < tol
            break;
        end
        x0 = x1;
        x1 = x2;
    end
    root = rad2deg(x2);
end

% ── Modified Secant Method ───────────────────────────────────────────────
function [root, hist, ea_hist, n] = mod_secant_method(f, x0d, delta, tol, maxiter)
% MOD_SECANT_METHOD  Modified secant method with single initial guess.
    x    = deg2rad(x0d);
    hist    = x0d;
    ea_hist = NaN;
    n       = 0;
    xn      = x;    % pre-assign so root is defined if maxiter=0

    for i = 1:maxiter
        n     = i;
        denom = f(x + delta*x) - f(x);
        if abs(denom) < eps
            warning('mod_secant_method: near-zero denominator at iteration %d. Stopping.', i);
            break;
        end
        xn = x - (delta * x * f(x)) / denom;
        ea = abs((xn - x) / xn) * 100;
        hist(end+1)    = rad2deg(xn);   %#ok<AGROW>
        ea_hist(end+1) = ea;            %#ok<AGROW>
        if ea < tol
            break;
        end
        x = xn;
    end
    root = rad2deg(xn);
end

% ── Newton-Raphson Method ────────────────────────────────────────────────
function [root, hist, ea_hist, n] = newton_raphson(f, fp, x0d, tol, maxiter)
% NEWTON_RAPHSON  Newton-Raphson method with analytical derivative.
    x    = deg2rad(x0d);
    hist    = x0d;
    ea_hist = NaN;
    n       = 0;
    xn      = x;    % pre-assign so root is defined if maxiter=0

    for i = 1:maxiter
        n = i;
        if abs(fp(x)) < eps
            error('newton_raphson: zero derivative at iteration %d. Choose a different initial guess.', i);
        end
        xn = x - f(x) / fp(x);
        ea = abs((xn - x) / xn) * 100;
        hist(end+1)    = rad2deg(xn);   %#ok<AGROW>
        ea_hist(end+1) = ea;            %#ok<AGROW>
        if ea < tol
            break;
        end
        x = xn;
    end
    root = rad2deg(xn);
end

% ── Trajectory Function ──────────────────────────────────────────────────
function y = traj(th_deg, x_arr, g, v0, y0)
% TRAJ  Compute ball trajectory y given launch angle (degrees) and x array.
    th = deg2rad(th_deg);
    y  = tan(th) .* x_arr ...
       - (g / (2 * v0^2 * cos(th)^2)) .* x_arr.^2 ...
       + y0;
end