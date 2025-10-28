function dx = coupledPendulums(~, x, p)
% coupledPendulums: RHS for a cart with two pendulums (planar).
% State (6×1): x = [phi1; dphi1; phi2; dphi2; q; dq]
%   phi1, phi2  : pendulum angles from vertical [rad]
%   q           : cart horizontal position [m]
%   d(·)        : time derivative
%
% Parameters in p (all optional; defaults below):
%   m1, m2  : bob masses [kg]
%   M       : cart mass [kg]
%   l       : pendulum length [m]
%   g       : gravity [m/s^2]
%   c_cart  : cart linear damping [N·s/m]
%   c_piv   : pivot viscous damping [N·m·s/rad] each pendulum

if nargin < 3 || isempty(p), p = struct(); end
if ~isfield(p,'m1'),     p.m1     = 0.1;  end
if ~isfield(p,'m2'),     p.m2     = 0.1;  end
if ~isfield(p,'M'),      p.M      = 5.0;  end
if ~isfield(p,'l'),      p.l      = 0.2;  end
if ~isfield(p,'g'),      p.g      = 9.81; end
if ~isfield(p,'c_cart'), p.c_cart = 0.2;  end
if ~isfield(p,'c_piv'),  p.c_piv  = 0.0;  end

% Unpack state
phi1 = x(1);   dphi1 = x(2);
phi2 = x(3);   dphi2 = x(4);
q    = x(5);   dq    = x(6);

% Trig shortcuts
s1 = sin(phi1);  c1 = cos(phi1);
s2 = sin(phi2);  c2 = cos(phi2);

% Convenience factors from the cart equation
alpha1 = p.m1*p.l/(p.M + p.m1 + p.m2);
alpha2 = p.m2*p.l/(p.M + p.m1 + p.m2);

% Denominator after eliminating ddphi1, ddphi2 via pendulum eqs
denom = 1 - c1^2*alpha1/p.l - c2^2*alpha2/p.l;
if abs(denom) < 1e-12, denom = sign(denom)*1e-12; end

% Numerator parts
num_core     = alpha1*( s1*dphi1^2 + (p.g/p.l)*c1*s1 ) + ...
               alpha2*( s2*dphi2^2 + (p.g/p.l)*c2*s2 );
num_cartDamp = -(p.c_cart/(p.M+p.m1+p.m2))*dq;
num_pivDamp  = (p.c_piv/(p.m1*p.l^2))*alpha1*c1*dphi1 + ...
               (p.c_piv/(p.m2*p.l^2))*alpha2*c2*dphi2;

ddq = (num_core + num_cartDamp + num_pivDamp) / denom;

% Angular accelerations:
ddphi1 = -(p.g/p.l)*s1 - (ddq/p.l)*c1 - (p.c_piv/(p.m1*p.l^2))*dphi1;
ddphi2 = -(p.g/p.l)*s2 - (ddq/p.l)*c2 - (p.c_piv/(p.m2*p.l^2))*dphi2;

% Assemble derivative
dx = [dphi1; ddphi1; dphi2; ddphi2; dq; ddq];
end
