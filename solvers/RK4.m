function [t, y] = RK4(fun, tspan, y0, h)
% RK4  Classic 4th-order Runge–Kutta for y' = f(t,y)
%   [t,y] = RK4(fun, [t0 tf], h, y0)
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

    k1 = fun(tk,       yk);
    k2 = fun(tk + h/2, yk + h*k1./2);
    k3 = fun(tk + h/2, yk + h*k2./2);
    k4 = fun(tk + h,   yk + h*k3);

    y(k+1,:) = (yk + (h/6)*(k1 + 2*k2 + 2*k3 + k4)).';

end
end