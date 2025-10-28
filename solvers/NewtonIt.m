function x = NewtonIt(G, x0, maxIt, tol)
% Solve G(x)=0 by Newton's method with numerical Jacobian.
% G   : function handle, R^n -> R^n
% x0  : initial guess
% maxIt (optional), tol (optional)

if nargin < 3 || isempty(maxIt), maxIt = 8;     end
if nargin < 4 || isempty(tol),   tol   = 1e-10; end

x = x0;
for it = 1:maxIt
    g  = G(x);
    J  = numjac(G, x);         % finite-difference Jacobian of G
    dx = - J \ g;              % solve J*dx = -g
    x  = x + dx;

    if norm(dx, inf) <= tol * (1 + norm(x, inf))
        break
    end
end
end

function J = numjac(F, x)
% Forward-difference Jacobian of vector function F at x
n  = numel(x);
J  = zeros(n);
Fx = F(x);
h  = sqrt(eps);     % ~1e-8 for double
for j = 1:n
    e      = zeros(n,1);
    e(j)   = h;
    J(:,j) = (F(x + e) - Fx) / h;
end
end
