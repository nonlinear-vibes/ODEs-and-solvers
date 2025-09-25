function dx = coupledPendulums(~, x)
% coupledPendulums:  RHS for a cart with two pendulums (planar).
% State vector (6×1):
%   x = [phi1; dphi1; phi2; dphi2; q; dq]
% where:
%   phi1, phi2  : pendulum angles from the vertical (radians)
%   q           : cart horizontal position
%   d(·)        : time derivative

% Dynamics parameters:
m1 = 0.1;     % pendulum bob mass
m2 = 0.1;     % pendulum bob mass
M  = 5;       % cart mass
l  = 0.2;     % pendulum length
g  = 9.81;    % gravity
c_cart = 0.2; % N·s/m   (cart linear damping)
c_piv  = 0;   % N·m·s/rad per pendulum (pivot damping)

% Unpack state
phi1 = x(1);   dphi1 = x(2);
phi2 = x(3);   dphi2 = x(4);
q    = x(5);   dq    = x(6);

% Trig shortcuts
s1 = sin(phi1);  c1 = cos(phi1);
s2 = sin(phi2);  c2 = cos(phi2);

% Convenience factor from the cart equation
%   ddq = α1 (s1*ω1^2 - c1*φ1'') + α2 (s2*ω2^2 - c2*φ2''),  
%   α1= m1*l/(M+m1+m2), α2= m2*l/(M+m1+m2)
alpha1 = m1*l/(M + m1 + m2);
alpha2 = m2*l/(M + m1 + m2);

% From pendulum equations:
%   φ1'' = -(g/l) s1 - (ddq/l) c1
%   φ2'' = -(g/l) s2 - (ddq/l) c2
% Substitute into the cart equation to get ddq:
denom = 1 - c1^2*alpha1/l - c2^2*alpha2/l;

% ddq = (alpha1 * (s1*dphi1^2 + (g/l)*c1*s1) + alpha2 * (s2*dphi2^2 + (g/l)*c2*s2) - (c_cart/(M+m1+m2))*dq) / denom;
% Numerator parts
num_core = alpha1*( s1*dphi1^2 + (g/l)*c1*s1 ) + ...
           alpha2*( s2*dphi2^2 + (g/l)*c2*s2 );

num_cartDamp = -(c_cart/(M+m1+m2))*dq;

% If pivot damping present, include its feedback into ddq:
num_pivDamp  = (c_piv/(m1*l^2))*alpha1*c1*dphi1 + ...
               (c_piv/(m2*l^2))*alpha2*c2*dphi2;

ddq = (num_core + num_cartDamp + num_pivDamp) / denom;

% Angular accelerations:
ddphi1 = -(g/l)*s1 - (ddq/l)*c1 - (c_piv/(m1*l^2))*dphi1;
ddphi2 = -(g/l)*s2 - (ddq/l)*c2 - (c_piv/(m2*l^2))*dphi2;

% Assemble derivative
dx = [dphi1; ddphi1; dphi2; ddphi2; dq; ddq];
end