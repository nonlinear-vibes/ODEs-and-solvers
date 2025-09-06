function [t, y] = EE(fun, tspan, y0, h)
% EE  Explicit Euler solver for y' = f(t,y)
%   [t,y] = EE(fun, [t0 tf], h, y0)
% - fun: function handle f(t,y) returning a column vector
% - y0 : initial state (row or column)

% time grid
t = (tspan(1):h:tspan(2)).';
N = numel(t);

% ensure column state inside, preallocate output as rows
m  = numel(y0);
y  = zeros(N, m);

y(1,:) = y0';           

% main loop
for k = 1:N-1
    tk = t(k);
    yk = y(k,:).';         % column state

    y(k+1,:) = yk + h*fun(tk, yk);
end
end

