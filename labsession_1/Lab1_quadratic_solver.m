% =========================================================
%  LAB 1 – Problem 1: Quadratic Equation Solver
%  Solves ax^2 + bx + c = 0 with full case handling.
%  Loops until user chooses to stop.
% =========================================================
clc; clear; close all;

while true
    % ----- User input -----
    a = input('Enter coefficient a: ');
    b = input('Enter coefficient b: ');
    c = input('Enter coefficient c: ');

    % ----- Special case: linear equation -----
    if a == 0
        if b == 0
            disp('Not an equation (a = b = 0).');
        else
            fprintf('Linear equation. Root: x = %.6f\n', -c/b);
        end
    else
        % ----- Discriminant -----
        disc = b^2 - 4*a*c;

        if disc > 0
            x1 = (-b + sqrt(disc)) / (2*a);
            x2 = (-b - sqrt(disc)) / (2*a);
            fprintf('Two real roots:\n  x1 = %.6f\n  x2 = %.6f\n', x1, x2);
        elseif disc == 0
            x1 = -b / (2*a);
            fprintf('One repeated real root: x = %.6f\n', x1);
        else
            realPart =  -b / (2*a);
            imagPart =  sqrt(-disc) / (2*a);
            fprintf('Complex roots:\n');
            fprintf('  x1 = %.6f + %.6fi\n', realPart,  imagPart);
            fprintf('  x2 = %.6f - %.6fi\n', realPart,  imagPart);
        end
    end

    % ----- Repeat? -----
    again = input('\nSolve another equation? (1 = Yes, 0 = No): ');
    if again ~= 1
        break;
    end
end
