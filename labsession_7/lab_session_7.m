clc; clear; close all;

if ~exist('graphs','dir'), mkdir('graphs'); end

%% ============================================================================
%  PROBLEM 1 -- Polynomial Differentiation
%% ============================================================================

%% 1A.  Function and exact analytical derivatives
% ─────────────────────────────────────────────────────────────────────────────
f  = @(x)  0.2 + 25*x  - 200*x.^2  + 675*x.^3  - 900*x.^4  + 400*x.^5;
f1 = @(x)  25  - 400*x + 2025*x.^2 - 3600*x.^3 + 2000*x.^4;
f2 = @(x) -400 + 4050*x - 10800*x.^2 + 8000*x.^3;
f3 = @(x)  4050 - 21600*x + 24000*x.^2;
f4 = @(x) -21600 + 48000*x;

x_exact = linspace(-1, 2, 1000);

%% 1B.  Plot f(x) and its four exact derivatives  (Figure 1)
% ─────────────────────────────────────────────────────────────────────────────
figure('Name','Exact Derivatives','NumberTitle','off');
funcs  = {f,  f1,  f2,  f3,  f4};
labels = {'f(x)', "f'(x)", "f''(x)", "f'''(x)", "f^{(4)}(x)"};
colors = {'b','r','g','m','k'};
for k = 1:5
    subplot(1,5,k);
    plot(x_exact, funcs{k}(x_exact), colors{k}, 'LineWidth', 2);
    yline(0,'k--','LineWidth',0.7);
    xlabel('x'); ylabel(labels{k}); title(labels{k});
    xlim([-1 2]); grid on;
end
sgtitle('f(x) and Its First Four Exact Derivatives', 'FontSize',11,'FontWeight','bold');
saveas(gcf,'graphs/lab7_fig1_exact_derivs.pdf');

%% 1C.  Numerical grid
% ─────────────────────────────────────────────────────────────────────────────
h  = 0.05;
x  = -1 : h : 2;          % same as linspace but exact spacing
n  = length(x);
fx = f(x);

%% 1D.  Forward-difference formulae  (loop-free, vectorised)
% Stencil applied pointwise; tail points filled with last valid value.
% Coefficients taken directly from Table 1 (image1 in the problem sheet).
% ─────────────────────────────────────────────────────────────────────────────

% -- 1st derivative ---------------------------------------------------------
d1_oh  = zeros(1,n);  d1_oh2 = zeros(1,n);
% O(h):  [-1  1] / h
d1_oh(1:n-1)   = (fx(2:n)   - fx(1:n-1))                           / h;
d1_oh(n)       = d1_oh(n-1);
% O(h^2): [-3  4  -1] / 2h
d1_oh2(1:n-2)  = (-fx(3:n)  + 4*fx(2:n-1) - 3*fx(1:n-2))          / (2*h);
d1_oh2(n-1:n)  = d1_oh2(n-2);

% -- 2nd derivative ---------------------------------------------------------
d2_oh  = zeros(1,n);  d2_oh2 = zeros(1,n);
% O(h):  [1  -2  1] / h^2
d2_oh(1:n-2)   = (fx(3:n)   - 2*fx(2:n-1) + fx(1:n-2))            / h^2;
d2_oh(n-1:n)   = d2_oh(n-2);
% O(h^2): [2  -5  4  -1] / h^2
d2_oh2(1:n-3)  = (2*fx(1:n-3) - 5*fx(2:n-2) + 4*fx(3:n-1) - fx(4:n)) / h^2;
d2_oh2(n-2:n)  = d2_oh2(n-3);

% -- 3rd derivative ---------------------------------------------------------
d3_oh  = zeros(1,n);  d3_oh2 = zeros(1,n);
% O(h):  [-1  3  -3  1] / h^3
d3_oh(1:n-3)   = (- fx(1:n-3) + 3*fx(2:n-2) - 3*fx(3:n-1) + fx(4:n)) / h^3;
d3_oh(n-2:n)   = d3_oh(n-3);
% O(h^2): [-5  18  -24  14  -3] / 2h^3
d3_oh2(1:n-4)  = (-5*fx(1:n-4) + 18*fx(2:n-3) - 24*fx(3:n-2) ...
                  + 14*fx(4:n-1) - 3*fx(5:n))                       / (2*h^3);
d3_oh2(n-3:n)  = d3_oh2(n-4);

% -- 4th derivative ---------------------------------------------------------
d4_oh  = zeros(1,n);  d4_oh2 = zeros(1,n);
% O(h):  [1  -4  6  -4  1] / h^4
d4_oh(1:n-4)   = (fx(1:n-4) - 4*fx(2:n-3) + 6*fx(3:n-2) ...
                  - 4*fx(4:n-1) + fx(5:n))                          / h^4;
d4_oh(n-3:n)   = d4_oh(n-4);
% O(h^2): [3  -14  26  -24  11  -2] / h^4
d4_oh2(1:n-5)  = (3*fx(1:n-5) - 14*fx(2:n-4) + 26*fx(3:n-3) ...
                  - 24*fx(4:n-2) + 11*fx(5:n-1) - 2*fx(6:n))       / h^4;
d4_oh2(n-4:n)  = d4_oh2(n-5);

%% 1E.  Print numerical comparison table at x = 0.5
% ─────────────────────────────────────────────────────────────────────────────
xs = 0.5;
[~, si] = min(abs(x - xs));
exact_v = [f1(xs), f2(xs), f3(xs), f4(xs)];
oh_v    = [d1_oh(si), d2_oh(si), d3_oh(si), d4_oh(si)];
oh2_v   = [d1_oh2(si), d2_oh2(si), d3_oh2(si), d4_oh2(si)];

fprintf('\n===== Numerical Derivatives at x = %.1f, h = %.2f =====\n', xs, h);
fprintf('%-12s  %12s  %12s  %12s  %12s  %12s\n', ...
        'Deriv', 'Exact', 'O(h)', 'Err O(h)', 'O(h^2)', 'Err O(h^2)');
fprintf('%s\n', repmat('-',1,75));
dnames = {'f1','f2','f3','f4'};
for k = 1:4
    fprintf('%-12s  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f\n', ...
            dnames{k}, exact_v(k), oh_v(k), abs(oh_v(k)-exact_v(k)), ...
            oh2_v(k), abs(oh2_v(k)-exact_v(k)));
end

%% 1F.  Figures 2–5: Exact vs O(h) vs O(h^2)
% ─────────────────────────────────────────────────────────────────────────────
fexact_cells = {f1,    f2,     f3,      f4     };
oh_cells     = {d1_oh, d2_oh,  d3_oh,   d4_oh  };
oh2_cells    = {d1_oh2,d2_oh2, d3_oh2,  d4_oh2 };
dtitles      = {'First Derivative',  'Second Derivative', ...
                'Third Derivative',  'Fourth Derivative'};
ylabels      = {"f'(x)", "f''(x)", "f'''(x)", "f^{(4)}(x)"};
fignums      = [2 3 4 5];
fignames     = {'lab7_fig2_deriv1','lab7_fig3_deriv2', ...
                'lab7_fig4_deriv3','lab7_fig5_deriv4'};

for k = 1:4
    figure('Name',sprintf('Deriv%d Num',k),'NumberTitle','off');
    plot(x_exact, fexact_cells{k}(x_exact), 'k-',  'LineWidth',2.4, ...
         'DisplayName','Exact'); hold on;
    plot(x,       oh_cells{k},              'r--', 'LineWidth',1.8, ...
         'DisplayName','O(h)');
    plot(x,       oh2_cells{k},             'b-.', 'LineWidth',1.8, ...
         'DisplayName','O(h^2)');
    yline(0,'Color',[.6 .6 .6],'LineWidth',0.6);
    xlabel('x'); ylabel(ylabels{k});
    title(sprintf('%s: Exact vs O(h) vs O(h^2)  (h = %.2f)', dtitles{k}, h));
    xlim([-1 2]); legend('Location','best'); grid on;
    saveas(gcf, sprintf('graphs/%s.pdf', fignames{k}));
end

%% 1G.  MATLAB built-ins: diff and gradient
% ─────────────────────────────────────────────────────────────────────────────
% diff(y)./diff(x) — O(h), length shrinks by 1 each call
% gradient(y, h)   — O(h^2) at interior, length preserved

% Pre-build diff and gradient results for all four orders
x_diff_cell = cell(1,4);  y_diff_cell = cell(1,4);
y_grad_cell = cell(1,4);

y_d = fx;  x_d = x;
y_g = fx;
for k = 1:4
    % diff
    y_d        = diff(y_d) ./ diff(x_d);
    x_d        = (x_d(1:end-1) + x_d(2:end)) / 2;   % midpoints
    x_diff_cell{k} = x_d;
    y_diff_cell{k} = y_d;
    % gradient
    y_g            = gradient(y_g, h);
    y_grad_cell{k} = y_g;
end

fignames_bi = {'lab7_fig6_diff_grad_deriv1','lab7_fig7_diff_grad_deriv2', ...
               'lab7_fig8_diff_grad_deriv3','lab7_fig9_diff_grad_deriv4'};
for k = 1:4
    figure('Name',sprintf('Deriv%d Builtin',k),'NumberTitle','off');
    plot(x_exact, fexact_cells{k}(x_exact), 'k-',  'LineWidth',2.4, ...
         'DisplayName','Exact'); hold on;
    plot(x_diff_cell{k}, y_diff_cell{k},    'r--', 'LineWidth',1.8, ...
         'DisplayName','diff');
    plot(x,              y_grad_cell{k},    'b-.', 'LineWidth',1.8, ...
         'DisplayName','gradient');
    yline(0,'Color',[.6 .6 .6],'LineWidth',0.6);
    xlabel('x'); ylabel(ylabels{k});
    title(sprintf('%s: Exact vs diff vs gradient', dtitles{k}));
    xlim([-1 2]); legend('Location','best'); grid on;
    saveas(gcf, sprintf('graphs/%s.pdf', fignames_bi{k}));
end

%% 1H.  Error vs step-size h at x = 1  (Figure 13)
% ─────────────────────────────────────────────────────────────────────────────
x_s      = 1.0;
exact_d1 = f1(x_s);
h_arr    = logspace(-8, 0, 200);

err_oh   = abs((f(x_s + h_arr) - f(x_s)) ./ h_arr - exact_d1);
err_oh2  = abs((-f(x_s+2*h_arr) + 4*f(x_s+h_arr) - 3*f(x_s)) ./ (2*h_arr) - exact_d1);

figure('Name','Error vs h','NumberTitle','off');
loglog(h_arr, err_oh,  'r-',  'LineWidth',2, 'DisplayName','O(h) formula'); hold on;
loglog(h_arr, err_oh2, 'b--', 'LineWidth',2, 'DisplayName','O(h^2) formula');
% Reference slopes
h_ref = [1e-5, 1e-1];
loglog(h_ref, 1e-8*(h_ref/1e-5).^1, 'k:',  'LineWidth',1.2, 'DisplayName','slope 1');
loglog(h_ref, 1e-10*(h_ref/1e-5).^2,'k--', 'LineWidth',1.2, 'DisplayName','slope 2');
xlabel('Step size h');  ylabel('Absolute error');
title('First-Derivative Error vs Step Size at x = 1');
legend('Location','northwest'); grid on;
saveas(gcf,'graphs/lab7_fig13_error_vs_h.pdf');

%% ============================================================================
%  PROBLEM 2 -- Radar Tracking  (Problem 21.14)
%% ============================================================================

%% 2A.  Data
% ─────────────────────────────────────────────────────────────────────────────
t_data  = [200, 202, 204, 206, 208, 210];
th_data = [0.75, 0.72, 0.70, 0.68, 0.67, 0.66];   % theta (rad)
r_data  = [5120, 5370, 5560, 5800, 6030, 6240];    % r (m)
h_t     = 2;       % time step (s)
idx     = 4;       % 1-based index for t = 206 s

%% 2B.  Centred O(h^2) differences at t = 206 s
% ─────────────────────────────────────────────────────────────────────────────
r_dot    = (r_data(idx+1)  - r_data(idx-1))  / (2*h_t);
r_ddot   = (r_data(idx+1)  - 2*r_data(idx)  + r_data(idx-1))  / h_t^2;
th_dot   = (th_data(idx+1) - th_data(idx-1)) / (2*h_t);
th_ddot  = (th_data(idx+1) - 2*th_data(idx) + th_data(idx-1)) / h_t^2;

r0 = r_data(idx);  th0 = th_data(idx);

%% 2C.  Velocity and acceleration in polar frame
% ─────────────────────────────────────────────────────────────────────────────
v_r   = r_dot;
v_th  = r0 * th_dot;
v_mag = sqrt(v_r^2 + v_th^2);

a_r   = r_ddot - r0 * th_dot^2;
a_th  = r0 * th_ddot + 2 * r_dot * th_dot;
a_mag = sqrt(a_r^2 + a_th^2);

fprintf('\n===== Problem 2 -- Centred Differences at t = 206 s =====\n');
fprintf('  r_dot         = %10.4f  m/s\n',    r_dot);
fprintf('  theta_dot     = %10.6f  rad/s\n',  th_dot);
fprintf('  r_ddot        = %10.4f  m/s^2\n',  r_ddot);
fprintf('  theta_ddot    = %10.6f  rad/s^2\n',th_ddot);
fprintf('\n  v_r           = %10.4f  m/s\n',   v_r);
fprintf('  v_theta       = %10.4f  m/s\n',    v_th);
fprintf('  |v|           = %10.4f  m/s\n',    v_mag);
fprintf('\n  a_r           = %10.4f  m/s^2\n', a_r);
fprintf('  a_theta       = %10.4f  m/s^2\n',  a_th);
fprintf('  |a|           = %10.4f  m/s^2\n',  a_mag);

%% 2D.  Cartesian conversion
% ─────────────────────────────────────────────────────────────────────────────
x_cart = r_data .* cos(th_data);
y_cart = r_data .* sin(th_data);

% Unit vectors at t = 206 s
er  = [cos(th0),  sin(th0)];
eth = [-sin(th0), cos(th0)];

% Velocity and acceleration in Cartesian at t = 206 s
v_cart = v_r * er  + v_th * eth;
a_cart = a_r * er  + a_th * eth;

x3 = x_cart(idx);  y3 = y_cart(idx);

%% 2E.  Figure 10 -- Rectangular path
% ─────────────────────────────────────────────────────────────────────────────
figure('Name','Radar Path','NumberTitle','off');
plot(x_cart, y_cart, 'bo-', 'LineWidth',2, 'MarkerFaceColor','b', ...
     'MarkerSize',7, 'DisplayName','Flight path');  hold on;
for k = 1:length(t_data)
    text(x_cart(k), y_cart(k), sprintf('  t=%d s', t_data(k)), ...
         'FontSize',8.5,'VerticalAlignment','bottom');
end
scatter(x3, y3, 120, 'r', 'pentagram', 'filled', ...
        'DisplayName','t = 206 s (analysis)');
xlabel('x (m)');  ylabel('y (m)');
title('Problem 2 -- Radar-Tracked Aircraft Path in Rectangular Coordinates');
legend('Location','northwest');  grid on;
saveas(gcf,'graphs/lab7_fig10_radar_path.pdf');

%% 2F.  Figure 11 -- Velocity vectors (at all interior points)
% ─────────────────────────────────────────────────────────────────────────────
scale_v = 3;   % display scale factor

figure('Name','Velocity Vectors','NumberTitle','off');
plot(x_cart, y_cart, 'bo-', 'LineWidth',1.8, 'MarkerFaceColor','b', ...
     'MarkerSize',5, 'Color',[0.3 0.5 0.9], 'DisplayName','Path');  hold on;

for k = 2:length(t_data)-1
    rd_k  = (r_data(k+1)  - r_data(k-1))  / (2*h_t);
    thd_k = (th_data(k+1) - th_data(k-1)) / (2*h_t);
    vr_k  = rd_k;
    vth_k = r_data(k) * thd_k;
    er_k  = [cos(th_data(k)),  sin(th_data(k))];
    eth_k = [-sin(th_data(k)), cos(th_data(k))];
    vc_k  = vr_k*er_k + vth_k*eth_k;
    quiver(x_cart(k), y_cart(k), vc_k(1)*scale_v, vc_k(2)*scale_v, 0, ...
           'r', 'LineWidth',2, 'MaxHeadSize',0.2, 'AutoScale','off', ...
           'DisplayName','');
end

scatter(x3, y3, 120, 'g', 'pentagram', 'filled', 'DisplayName','t = 206 s');
xlabel('x (m)');  ylabel('y (m)');
title('Problem 2 -- Velocity Vectors on the Aircraft Path');
legend({'Path','Velocity \vec{v}','t = 206 s'}, 'Location','northwest');
grid on;
saveas(gcf,'graphs/lab7_fig11_velocity.pdf');

%% 2G.  Figure 12 -- Acceleration vectors (at all interior points)
% ─────────────────────────────────────────────────────────────────────────────
scale_a = 25;

figure('Name','Acceleration Vectors','NumberTitle','off');
plot(x_cart, y_cart, 'bo-', 'LineWidth',1.8, 'MarkerFaceColor','b', ...
     'MarkerSize',5, 'Color',[0.3 0.5 0.9], 'DisplayName','Path');  hold on;

for k = 2:length(t_data)-1
    rd_k   = (r_data(k+1)  - r_data(k-1))  / (2*h_t);
    thd_k  = (th_data(k+1) - th_data(k-1)) / (2*h_t);
    rdd_k  = (r_data(k+1)  - 2*r_data(k)  + r_data(k-1))  / h_t^2;
    thdd_k = (th_data(k+1) - 2*th_data(k) + th_data(k-1)) / h_t^2;
    ar_k   = rdd_k  - r_data(k)*thd_k^2;
    ath_k  = r_data(k)*thdd_k + 2*rd_k*thd_k;
    er_k   = [cos(th_data(k)),  sin(th_data(k))];
    eth_k  = [-sin(th_data(k)), cos(th_data(k))];
    ac_k   = ar_k*er_k + ath_k*eth_k;
    quiver(x_cart(k), y_cart(k), ac_k(1)*scale_a, ac_k(2)*scale_a, 0, ...
           'Color',[0.9 0.4 0], 'LineWidth',2, 'MaxHeadSize',0.2, ...
           'AutoScale','off', 'DisplayName','');
end

scatter(x3, y3, 120, 'g', 'pentagram', 'filled', 'DisplayName','t = 206 s');
xlabel('x (m)');  ylabel('y (m)');
title('Problem 2 -- Acceleration Vectors on the Aircraft Path');
legend({'Path','Acceleration \vec{a}','t = 206 s'}, 'Location','northwest');
grid on;
saveas(gcf,'graphs/lab7_fig12_acceleration.pdf');

fprintf('\nAll figures saved to graphs/ folder.\n');
