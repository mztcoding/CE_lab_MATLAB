clc; clear; close all;

f = @(t, y) -0.1*exp(-0.1*t).*cos(2*t) - 2*exp(-0.1*t).*sin(2*t);
t_ic = 10; y_ic = 10;
t_0 = -5; t_end = 25;

h = input('Step size: ');

[t1, y1] = rk_method(f, t_ic, y_ic, t_0, t_end, h);

plot(t1, y1, 'b', 'LineWidth', 2);
xlabel('t'); ylabel('y(t)'); grid on;

function [t_all, y_all] = rk_method(f, t_ic, y_ic, t_0, t_end, h)
    % Forward integration
    t = (t_ic:h:t_end)';
    y = zeros(length(t), 1);
    y(1) = y_ic;
    for i = 1:length(t)-1
        k1 = f(t(i), y(i));
        k2 = f(t(i)+h/2, y(i)+h/2*k1);
        k3 = f(t(i)+h/2, y(i)+h/2*k2);
        k4 = f(t(i)+h, y(i)+h*k3);
        avg_slope = (k1 + 2*k2 + 2*k3 + k4) / 6;
        y(i+1) = y(i) + avg_slope * h;
    end

    % Backward integration
    t_back = (t_ic:-h:t_0)';
    y_back = zeros(length(t_back), 1);
    y_back(1) = y_ic;
    for i = 1:length(t_back)-1
        k1 = f(t_back(i), y_back(i));
        k2 = f(t_back(i)-h/2, y_back(i)-h/2*k1);
        k3 = f(t_back(i)-h/2, y_back(i)-h/2*k2);
        k4 = f(t_back(i)-h, y_back(i)-h*k3);
        avg_slope = (k1 + 2*k2 + 2*k3 + k4) / 6;
        y_back(i+1) = y_back(i) - avg_slope * h;
    end

    t_all = [flipud(t_back); t(2:end)];
    y_all = [flipud(y_back); y(2:end)];
end