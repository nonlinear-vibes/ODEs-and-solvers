function x = NewtonIt(fun, tkp1, x, yk, fk, h, maxIt)
% Newton iteration for G(x) = x - yk - (h/2)*( f(tk,yk) + f(tk+h,x) ) = 0
% with numerical Jacobian: dG = I - (h/2)*J(tk+h,x)
%
% Inputs:
%   fun  : RHS f(t,y)
%   tkp1 : time t_{i+1}
%   x    : initial guess for y_{i+1}
%   yk   : previous state y_k
%   fk   : f(t_i, y_i)
%   h    : step size
%   maxIt: max Newton iterations (no line search here)
%
% Output:
%   x    : Newton iterate (approx y_{i+1})

tol = 1e-10;          

for it = 1:maxIt

    fx = fun(tkp1, x);
    Jx = numjac(@(z) fun(tkp1, z), x);

    % residual and Newton matrix
    G  = x - yk - 0.5*h*(fk + fx);
    A  = eye(numel(x)) - 0.5*h*Jx;

    % Newton step
    dx = -A \ G;
    x  = x + dx;

    % Convergence check on the step
    if norm(dx, inf) <= tol * (1 + norm(x, inf))
        break
    end
end

end

function J = numjac(fun_y, x)
% Simple forward-difference Jacobian w.r.t. y
n  = numel(x);
J  = zeros(n);
fx = fun_y(x);
h  = sqrt(eps);
for j = 1:n
    e    = zeros(n,1);
    e(j) = h;
    J(:,j) = (fun_y(x + e) - fx) / h;
end
end