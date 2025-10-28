function [t, y] = TRAP(fun, tspan, y0, h)
% TRAP  Implicit trapezoidal method for y' = f(t,y)
% Usage: [t, y] = TRAP(@LV, tspan, h, y0)
% Solves the nonlinear step with Newton

% grid
t = (tspan(1):h:tspan(2)).';
N = numel(t);

% store states row-wise; pass columns to fun
m = numel(y0);
y = zeros(N, m);

y(1,:) = y0';

for k = 1:N-1
    tk   = t(k);
    tkp1 = t(k+1);
    yk   = y(k,:)';     % column vector for f

    % evaluate f at (tk, yk)
    fk = fun(tk, yk);
    
    % explicit predictor (Euler) for the intial guess of Newton
    x  = yk + h*fk;

    % Newton iterations at (t_{k+1}, x)
    G = @(x) x - yk - 0.5*h*(fk + fun(tkp1, x));
    x = NewtonIt(G, x);

    % trapezoid update (re-uses fk, evaluates f at the converged x)
    y(k+1,:) = (yk + h/2*(fk+fun(tkp1, x)))';
end
end
