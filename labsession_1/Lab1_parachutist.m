% =========================================================
%  LAB 1 – Problem 2: Falling Parachutist
%  Compares the analytical solution to Euler numerical
%  solutions for linear and second-order drag models.
%  m = 68.1 kg, c = 12.5 kg/s, c' = 0.225 kg/m
% =========================================================
clc; clear; close all;

%% ===== PARAMETERS =====
m   = 68.1;    % mass (kg)
c   = 12.5;    % linear drag coefficient (kg/s)
g   = 9.8;     % gravity (m/s^2)
dt  = 2;       % time step (s)
t_end = 14;    % end time (s)

%% ===== ANALYTICAL SOLUTION =====
t_anal = 0:0.1:t_end;
v_anal = (g*m/c) .* (1 - exp(-(c/m).*t_anal));

%% ===== NUMERICAL (EULER) -- LINEAR DRAG =====
t_num = 0:dt:t_end;
v_num = zeros(size(t_num));
% v(0) = 0 (starts from rest)
for i = 1:length(t_num)-1
    v_num(i+1) = v_num(i) + (g - (c/m)*v_num(i)) * dt;
end

fprintf('Euler (linear drag) results:\n');
fprintf('  t(s)    v_analytical    v_numerical\n');
for i = 1:length(t_num)
    v_a = (g*m/c)*(1 - exp(-(c/m)*t_num(i)));
    fprintf('  %4.1f    %10.4f    %10.4f\n', t_num(i), v_a, v_num(i));
end

%% ===== NUMERICAL (EULER) -- SECOND-ORDER DRAG (Prob 1.3) =====
c2    = 0.225;   % second-order drag coefficient (kg/m)
v_num2 = zeros(size(t_num));
for i = 1:length(t_num)-1
    v_num2(i+1) = v_num2(i) + (g - (c2/m)*v_num2(i)^2) * dt;
end

%% ===== PLOT =====
figure;
plot(t_anal, v_anal, 'k-',  'LineWidth', 2, 'DisplayName', 'Analytical (linear drag)');
hold on;
plot(t_num,  v_num,  'bo-', 'LineWidth', 1.5, 'MarkerFaceColor','b', ...
     'DisplayName', 'Euler numerical (linear drag)');
plot(t_num,  v_num2, 'rs-', 'LineWidth', 1.5, 'MarkerFaceColor','r', ...
     'DisplayName', 'Euler numerical (2nd-order drag)');
yline((g*m/c), 'k--', 'LineWidth', 1, 'DisplayName', 'Terminal velocity');
xlabel('Time, t (s)');
ylabel('Velocity, v (m/s)');
title('Falling Parachutist: Analytical vs. Numerical Solution');
legend('Location', 'southeast');
grid on;
