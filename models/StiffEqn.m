function y_prime = StiffEqn(t, y)
% STIFFEQN Nonlinear 1D ODE with stiff behaviour(!)
%   y' = lambda * ( t^alpha * y^(alpha+1) - 1/t ) - 1/t^2
%   The RHS is singular at t = 0 → start from t0 > 0 (e.g., t0 = 0.5).
%   Analytical solution: y(t) = 1/t (t>0) 

lambda = 0.8;
alpha  = 2;

y_prime = lambda* ((t^alpha)*y^(alpha+1)-1/t) - 1/t^2;

end

