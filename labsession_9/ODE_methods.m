clc; clear; close all;

h = input('Enter step size h: ');

t0 = 0; tf = 50; y0 = 1;
f  = @(t,y) y * sin(t);        % dy/dt = y*sin(t)

% Analytical Solution: y = e^(1 - cos(t))
ta = t0:0.01:tf;                % Fine grid for smooth curve
ya = exp(1 - cos(ta));          % Direct formula — exact solution

% 2. Euler's Method
t = t0:h:tf;  N = length(t);
y_euler = zeros(1,N);  y_euler(1) = y0;
for i = 1:N-1
    y_euler(i+1) = y_euler(i) + h * f(t(i), y_euler(i));
end

% 3. Midpoint Method
y_mid = zeros(1,N);  y_mid(1) = y0;
for i = 1:N-1
    k1 = f(t(i),       y_mid(i));
    k2 = f(t(i)+h/2,   y_mid(i)+h/2*k1);
    y_mid(i+1) = y_mid(i) + h*k2;
end

% 4. Predictor-Corrector
y_pc = zeros(1,N);  y_pc(1) = y0;
for i = 1:N-1
    yp           = y_pc(i) + h * f(t(i), y_pc(i));
    y_pc(i+1)    = y_pc(i) + h/2*(f(t(i),y_pc(i)) + f(t(i+1),yp));
end

% 5. Iterative Heun's Method
y_heun = zeros(1,N);  y_heun(1) = y0;
for i = 1:N-1
    f0 = f(t(i), y_heun(i));
    yp = y_heun(i) + h*f0;
    for k = 1:50
        yp_new = y_heun(i) + h/2*(f0 + f(t(i+1),yp));
        if abs(yp_new - yp) < 1e-6, break; end
        yp = yp_new;
    end
    y_heun(i+1) = yp_new;
end

% Plot
figure; hold on;
plot(ta, ya,      'k-',  'LineWidth', 2,   'DisplayName', 'Analytical');
plot(t,  y_euler, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Euler');
plot(t,  y_mid,   'b:',  'LineWidth', 2,   'DisplayName', 'Midpoint');
plot(t,  y_pc,    'g-.', 'LineWidth', 1.5, 'DisplayName', 'Predictor-Corrector');
plot(t,  y_heun,  'm-',  'LineWidth', 1.5, 'DisplayName', 'Iterative Heun''s');
xlabel('t'); ylabel('y(t)');
title(sprintf('ODE Methods Comparison  (h = %g)', h));
legend; grid on;