function f = VDP(~, y)
% VDP  van der Pol oscillator with inline parameter mu.
% State y = [y1; y2].  y1' = y2, y2' = mu*(1 - y1^2)*y2 - y1

% Nonlinear damping coeff.
mu = 3;
a  = 1; 
b  = 1;

% System
f  = [ y(2);
       mu*(a - y(1)^2)*y(2) - b*y(1) ];

end

