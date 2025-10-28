function y_prime = StiffEqn(t, y, p)
% STIFFEQN  Nonlinear 1D ODE with stiff behavior.
%
%   y' = λ( t^α * y^(α+1) - 1/t ) - 1/t^2
%
%   The RHS is singular at t = 0 (1/t and 1/t^2 terms), so start from t0 > 0,
%   e.g., t0 = 0.5 or 1.0.
%
%   Analytical solution (for λ = α + 1):  y(t) = 1/t for t > 0.
%
% Parameters in p (optional):
%   p.lambda   (default 0.8)
%   p.alpha    (default 2)
%
% Usage:
%   p = struct('lambda', 0.5, 'alpha', 3);
%   f = @(t,y) StiffEqn(t,y,p);
%   [t, y] = IE(f, [0.5, 10], 2, 0.01);

% Handle optional parameter struct
if nargin < 3 || isempty(p), p = struct(); end
if ~isfield(p, 'lambda'), p.lambda = 0.8; end
if ~isfield(p, 'alpha'),  p.alpha  = 2;   end

λ = p.lambda;
α = p.alpha;

% RHS
y_prime = λ*( (t^α)*y^(α+1) - 1/t ) - 1/t^2;

end
