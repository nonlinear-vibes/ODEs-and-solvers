function f = LV(~, y)
% LV  Lotka–Volterra predator–prey system.
% State y = [x; z] where x = prey, z = predator. Returns a column vector.
%
%   x' = a x - b x z
%   z' = c x z - d z

% System parameters (edit these)
a = 2/3;    % prey growth
b = 4/3;    % predation rate
c = 1;      % predator growth per prey eaten
d = 1;      % predator death

% System
f = [ a*y(1) - b*y(1)*y(2);
      c*y(1)*y(2) - d*y(2) ];

end

