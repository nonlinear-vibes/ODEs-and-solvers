function dx = CR3BP(~, x, p)
% CR3BP: Circular restricted three-body problem in rotating frame (ω=1).
% State x = [X; VX; Y; VY]. Primaries at (-mu,0) and ((1-mu),0).
%
% Params:
%   p.mu ∈ (0, 0.5)  mass parameter (default Earth–Moon: 0.012277471)
%
% Equations:
%   X'  = VX
%   Y'  = VY
%   VX' =  X + 2*VY - nu*(X+mu)/r1^3 - mu*(X-nu)/r2^3
%   VY' =  Y - 2*VX - nu*Y/r1^3     - mu*Y/r2^3
%
% Here we use the Earth–Moon mass ratio:
%   mu = 0.012277471 (Moon), nu = 1 - mu (Earth).
%
% Inputs:
%   ~ : time (unused; system is autonomous in the rotating frame)
%   x : 4x1 state [X; VX; Y; VY]
% Output:
%   dx: 4x1 time derivative
%
% Three known periodic orbits:
% 1) Tspan = [0, 6.2];   x0 = [1.2; 0; 0; -1.049357510];
% 2) Tspan = [0, 11.2];  x0 = [0.994; 0; 0; -2.0317326295573368];
% 3) Tspan = [0, 17.1];  x0 = [0.994; 0; 0; -2.0015851063790825];

if nargin < 3 || isempty(p), p = struct(); end
if ~isfield(p,'mu'), p.mu = 0.012277471; end
nu = 1 - p.mu;

% Unpack state
X  = x(1);  VX = x(2);
Y  = x(3);  VY = x(4);

% Distances to the primaries
r1sq = (X + p.mu)^2 + Y^2;      % distance^2 to (-mu, 0)
r2sq = (X - nu)^2 + Y^2;      % distance^2 to (+nu, 0)
r1_3 = r1sq^(3/2);
r2_3 = r2sq^(3/2);

% Dynamics in rotating frame
aX = X + 2*VY - nu*(X + p.mu)/r1_3 - p.mu*(X - nu)/r2_3;
aY = Y - 2*VX - nu*Y/r1_3       - p.mu*Y/r2_3;

% Derivative
dx = [VX; aX; VY; aY];
end
