function [xr, roots_arr, ea_arr, n_iter] = false_position(f, xl, xu, tol, max_iter)
% FALSE_POSITION  Regula Falsi root-finding method.
%
%   USAGE:
%     [xr, roots_arr, ea_arr, n_iter] = false_position(f, xl, xu, tol, max_iter)
%
%   INPUTS:
%     f        - anonymous function handle, e.g. @(x) x^2 - 4
%     xl, xu   - initial bracket (must satisfy f(xl)*f(xu) < 0)
%     tol      - stopping tolerance (approximate % relative error)
%     max_iter - maximum number of iterations
%
%   OUTPUTS:
%     xr        - root estimate
%     roots_arr - array of all root estimates per iteration
%     ea_arr    - array of approximate % relative errors per iteration
%     n_iter    - total iterations performed
%
%   NOTE: Save this file as false_position.m in the same folder as
%         Lab3_fp_x10.m, Lab3_fp_rocket.m, and Lab3_fp_sphere.m

if f(xl)*f(xu) >= 0
    error('false_position: f(xl)*f(xu) >= 0. No sign change detected.');
end

xr_old    = xl;
ea        = inf;
roots_arr = [];
ea_arr    = [];
n_iter    = 0;

for i = 1:max_iter
    n_iter = i;
    xr     = xu - f(xu)*(xl - xu) / (f(xl) - f(xu));

    if xr ~= 0
        ea = abs((xr - xr_old) / xr) * 100;
    end

    roots_arr(end+1) = xr;
    ea_arr(end+1)    = ea;

    if f(xl)*f(xr) < 0
        xu = xr;
    elseif f(xr)*f(xu) < 0
        xl = xr;
    else
        break;
    end

    if ea < tol,  break;  end
    xr_old = xr;
end
end
