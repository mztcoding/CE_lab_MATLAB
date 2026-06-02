
%  NUMERICAL INTEGRATION — Problems 19.8a and 19.6

clc; clear; close all;


%  PROBLEM 19.8a — Distance from Velocity Data (Trapezoidal Rule)

fprintf('--- PROBLEM 19.8a: Distance from Velocity Data ---\n\n');

t_data = [1, 2, 3.25, 4.5, 6, 7, 8, 8.5, 9, 10];
v_data = [5, 6, 5.5,  7,  8.5, 8, 6, 7,  7,  5];

distance_trap = trapz(t_data, v_data);
T_total       = t_data(end) - t_data(1);
avg_velocity  = distance_trap / T_total;

fprintf('Total Distance Traveled = %.4f units\n', distance_trap);
fprintf('Total Time Interval     = %.4f units\n', T_total);
fprintf('Average Velocity        = %.4f units/time\n\n', avg_velocity);

% --- Plot 19.8a ---
figure('Name','Problem 19.8a - Velocity Data','NumberTitle','off',...
    'Color','w','Position',[100 100 850 460]);

hold on; grid on; box on;
set(gca,'Color','w','GridColor',[0.85 0.85 0.85],'GridAlpha',1,'FontSize',11);

for i = 1:length(t_data)-1
    fill([t_data(i), t_data(i+1), t_data(i+1), t_data(i)], ...
         [0, 0, v_data(i+1), v_data(i)], ...
         [0.75 0.90 1.0], 'EdgeColor',[0.5 0.5 0.5], ...
         'LineWidth', 0.8, 'FaceAlpha', 0.6, 'HandleVisibility','off');
end

plot(t_data, v_data, 'b-o', 'LineWidth', 2, ...
    'MarkerFaceColor','b', 'MarkerSize', 8, 'DisplayName','Velocity Data');

yline(avg_velocity, 'r--', 'LineWidth', 1.8, 'DisplayName', ...
    sprintf('Avg Velocity = %.4f', avg_velocity));

for i = 1:length(t_data)-1
    seg_area = (t_data(i+1) - t_data(i)) * (v_data(i) + v_data(i+1)) / 2;
    xm = (t_data(i) + t_data(i+1)) / 2;
    ym = (v_data(i) + v_data(i+1)) / 4;
    text(xm, ym, sprintf('%.2f', seg_area), 'FontSize', 8.5, ...
        'HorizontalAlignment','center','Color',[0.15 0.15 0.55],'FontWeight','bold');
end

xlabel('Time t', 'FontSize', 12, 'FontWeight','bold');
ylabel('Velocity v', 'FontSize', 12, 'FontWeight','bold');
title(sprintf('Problem 19.8a — Trapezoidal Rule on Velocity Data\nTotal Distance = %.4f  |  Avg Velocity = %.4f', ...
    distance_trap, avg_velocity), 'FontSize', 11, 'FontWeight','bold');
legend('Location','northwest','FontSize',10);
xlim([0.5, 10.5]); ylim([0, 11]);



%  PROBLEM 19.6 — Double Integral (Composite Trapezoidal Rule)
%  Integral: x in [-2,2], y in [0,4] of (x^2 - 3y^2 + x*y^3) dx dy

fprintf('--- PROBLEM 19.6: Double Integral ---\n');
fprintf('f(x,y) = x^2 - 3y^2 + x*y^3\n');
fprintf('Limits: x in [-2, 2],  y in [0, 4]\n\n');

n = input('Enter n for composite trapezoidal rule (must be even, default = 2): ');
if isempty(n) || mod(n, 2) ~= 0
    fprintf('Using default n = 2\n');
    n = 2;
end

% Correct limits: inner integral over x in [0,4], outer over y in [-2,2]
x1 =  0;  x2 = 4;
y1 = -2;  y2 = 2;
f  = @(x, y) x.^2 - 3*y.^2 + x.*y.^3;

exact_19_6 = 64/3;

% Integrate in x first (inner), then y (outer)
y_nodes = linspace(y1, y2, n+1);
I_x = zeros(1, n+1);
for j = 1:n+1
    I_x(j) = composite_trapezoid(@(x) f(x, y_nodes(j)), x1, x2, n);
end
I_trap   = composite_trapezoid_vec(I_x, y1, y2, n);
eps_trap = abs((exact_19_6 - I_trap) / exact_19_6) * 100;

fprintf('\nResults for Problem 19.6:\n');
fprintf('  Analytical Result              = %.6f\n', exact_19_6);
fprintf('  Composite Trapezoidal (n = %d) = %.6f  |  eps_t = %.4f%%\n', ...
    n, I_trap, eps_trap);

% --- Plot 19.6 ---
figure('Name','Problem 19.6 - Double Integral','NumberTitle','off',...
    'Color','w','Position',[120 120 680 520]);

hold on; grid on; box on;
set(gca,'Color','w','GridColor',[0.85 0.85 0.85],'GridAlpha',1,'FontSize',11);

x_surf = linspace(x1, x2, 50);
y_surf = linspace(y1, y2, 50);
[X_s, Y_s] = meshgrid(x_surf, y_surf);
Z_s = f(X_s, Y_s);
surf(X_s, Y_s, Z_s, 'EdgeAlpha', 0.15, 'FaceAlpha', 0.85);
colormap(parula); colorbar;

x_nodes = linspace(x1, x2, n+1);
y_nodes_plot = linspace(y1, y2, n+1);
[Xn, Yn] = meshgrid(x_nodes, y_nodes_plot);
Zn = f(Xn, Yn);
plot3(Xn(:), Yn(:), Zn(:), 'ro', 'MarkerFaceColor','r', ...
    'MarkerSize', 7, 'DisplayName', sprintf('Trap Nodes (n=%d)', n));

xlabel('x  [0, 4]', 'FontSize', 11, 'FontWeight','bold');
ylabel('y  [-2, 2]', 'FontSize', 11, 'FontWeight','bold');
zlabel('f(x,y)', 'FontSize', 11, 'FontWeight','bold');
title(sprintf('Problem 19.6 — f(x,y) = x^2 - 3y^2 + xy^3\nComposite Trap (n=%d): %.6f  |  Exact: %.6f  |  \\epsilon_t = %.4f%%', ...
    n, I_trap, exact_19_6, eps_trap), 'FontSize', 10, 'FontWeight','bold');
legend('Location','best','FontSize', 9);
view([-35 30]);



%  LOCAL FUNCTIONS


function I = composite_trapezoid(f, a, b, n)
    x = linspace(a, b, n+1);
    y = f(x);
    h = (b - a) / n;
    I = h * (y(1)/2 + sum(y(2:end-1)) + y(end)/2);
end

function I = composite_trapezoid_vec(y, a, b, n)
    h = (b - a) / n;
    I = h * (y(1)/2 + sum(y(2:end-1)) + y(end)/2);
end