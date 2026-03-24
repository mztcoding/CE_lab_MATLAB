% =========================================================
%  LAB 4 – Fixed-Point Iteration
%  Spring-mass deflection: k1=40000 N/m, k2=40 N/m^1.5
%  m=95 kg, g=9.81 m/s^2, h=0.43 m
%  g(d) = sqrt((mg*d + mgh - (2k2/5)*d^2.5) / (0.5*k1))
%  Enter initial guess when prompted (try 0.1).
% =========================================================
clc; clear;

%% Parameters
k1=40000; k2=40; m=95; g=9.81; h=0.43;
mg=m*g; mgh=mg*h;
g_fun = @(d) sqrt((mg.*d + mgh - (2*k2/5).*d.^(2.5))./(0.5*k1));

%% Iteration
d = input('Enter initial guess (try 0.1): ');
tol=0.001; errors=[];
for i=1:100
    d_new=g_fun(d); errors(end+1)=abs(d_new-d);
    fprintf('Iter %d: d = %.8f\n',i,d_new);
    if errors(end)<tol, break; end
    d=d_new;
end
root=d_new;
fprintf('Root = %.8f m\n',root);

%% Plot 1 – convergence
figure; plot(1:length(errors),errors,'-o');
xlabel('Iteration'); ylabel('|d_{new}-d|');
title('Fixed-Point Convergence'); grid on;

%% Plot 2 – g(d) vs y = d (cobweb), range covers root ~0.167
d_vals = linspace(0, 0.25, 500);
figure; plot(d_vals,g_fun(d_vals),'r','LineWidth',1.5); hold on;
plot(d_vals,d_vals,'b--','LineWidth',1.2);
plot(root,root,'ko','MarkerFaceColor','y','MarkerSize',10);
text(root+0.005,root,sprintf('Fixed Point\n(%.5f)',root),'FontSize',9);
xlabel('d (m)'); ylabel('Value');
legend('g(d)','y = d','Fixed Point','Location','northwest');
title('Graphical Fixed-Point Method – Spring Mass'); grid on;
