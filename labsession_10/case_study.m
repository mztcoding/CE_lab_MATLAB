clc; clear; close all;

h = 0.0625; tspan = [0 40]; y0 = [2 1];
a = 1.2; b = 0.6; c = 0.8; d = 0.3;

[t, y] = eulersys(@predprey, tspan, y0, h, a, b, c, d);
subplot(2,2,1); plot(t, y(:,1), t, y(:,2), '--')
legend('prey','predator'); title('(a) Euler time plot')
subplot(2,2,2); plot(y(:,1), y(:,2))
title('(b) Euler phase plane plot')

[t, y] = rk4sys(@predprey, tspan, y0, h, a, b, c, d);
subplot(2,2,3); plot(t, y(:,1), t, y(:,2), '--')
title('(c) RK4 time plot')
subplot(2,2,4); plot(y(:,1), y(:,2))
title('(d) RK4 phase plane plot')

h = 0.03125; tspan = [0 20]; y0 = [5 5 5];
sigma = 10; b = 8/3; r = 28;

[t, y] = rk4sys(@lorenz, tspan, y0, h, sigma, b, r);
figure;
plot(t, y(:,1)); title('Lorenz x vs t')
xlabel('t'); ylabel('x'); grid on

figure;
subplot(1,3,1); plot(y(:,1), y(:,2)); xlabel('x'); ylabel('y'); title('(a) y vs x')
subplot(1,3,2); plot(y(:,1), y(:,3)); xlabel('x'); ylabel('z'); title('(b) z vs x')
subplot(1,3,3); plot(y(:,2), y(:,3)); xlabel('y'); ylabel('z'); title('(c) z vs y')

figure;
plot3(y(:,1), y(:,2), y(:,3))
xlabel('x'); ylabel('y'); zlabel('z'); grid on

function yp = predprey(t, y, a, b, c, d)
    yp = [a*y(1) - b*y(1)*y(2); -c*y(2) + d*y(1)*y(2)];
end

function yp = lorenz(t, y, sigma, b, r)
    yp = [-sigma*y(1) + sigma*y(2); r*y(1) - y(2) - y(1)*y(3); -b*y(3) + y(1)*y(2)];
end

function [t, y] = eulersys(f, tspan, y0, h, varargin)
    t = (tspan(1):h:tspan(2))';
    y = zeros(length(t), length(y0));
    y(1,:) = y0;
    for i = 1:length(t)-1
        y(i+1,:) = y(i,:) + h * f(t(i), y(i,:)', varargin{:})';
    end
end

function [t, y] = rk4sys(f, tspan, y0, h, varargin)
    t = (tspan(1):h:tspan(2))';
    y = zeros(length(t), length(y0));
    y(1,:) = y0;
    for i = 1:length(t)-1
        k1 = f(t(i), y(i,:)', varargin{:});
        k2 = f(t(i)+h/2, y(i,:)'+h/2*k1, varargin{:});
        k3 = f(t(i)+h/2, y(i,:)'+h/2*k2, varargin{:});
        k4 = f(t(i)+h, y(i,:)'+h*k3, varargin{:});
        avg_slope = (k1 + 2*k2 + 2*k3 + k4) / 6;
        y(i+1,:) = y(i,:) + (avg_slope * h)';
    end
end