function [t, y] = IE(fun, tspan, y0, h) 
% IE  Implicit Euler (backward Euler) for y' = f(t,y)
% Usage: [t, y] = IE(@LV, tspan, h, y0)
% Solves the nonlinear step with Newton

% grid
t = (tspan(1):h:tspan(2)).';
N = numel(t);

% store states row-wise; pass columns to fun
m  = numel(y0);
y  = zeros(N, m);
y(1,:) = y0';

for k = 1:N-1
    tk   = t(k);
    tkp1 = t(k+1);
    yk   = y(k,:)';

    fk = fun(tk, yk);
    x  = yk + h*fk;     % Euler predictor

    % Solve G_BE(x) = x - yk - h f(t_{k+1}, x) = 0 via NewtonIt
    G = @(x) x - yk - h * fun(tkp1, x);
    x = NewtonIt(G, x);

    % implicit Euler update is just the converged x
    % transpose to store as a row
    y(k+1,:) = x.';     
end
end



