function f = Rayleigh(~,y)
% RAYLEIGH  Rayleigh oscillator (nonlinear damping on y')
%   y1' = y2
%   y2' = mu * (1 - y2^2) * y2 - y1
%
% This is the “Rayleigh equation” form: y'' + y = mu * (1 - (y')^2) * y'
% (contrast with van der Pol: y'' + y = mu*(1 - y^2)*y')

% nonlinearity / damping strength
mu = 3;

% RHS
f = [ y(2);
      mu*(1 - y(2)^2)*y(2) - y(1) ];

end