function f = SIR(~, y, p)
% SIR: Basic SIR compartmental epidemic model.
%   S' = -beta * S * I
%   I' =  beta * S * I - gamma * I
%   R' =  gamma * I

if nargin < 3 || isempty(p), p = struct(); end
if ~isfield(p,'beta'),  p.beta  = 0.30; end  % infection rate
if ~isfield(p,'gamma'), p.gamma = 0.10; end  % recovery rate

S = y(1); I = y(2); R = y(3);
lambda = p.beta * S * I;

f = [ -lambda;
       lambda - p.gamma * I;
       p.gamma * I ];
end

