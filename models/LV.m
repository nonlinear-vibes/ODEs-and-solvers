function f = LV(~, y)
% LV: Lotka–Volterra predator–prey system.
% y = [x; z] where x = prey, z = predator.
%   x' = a x - b x z
%   z' = c x z - d z

if nargin < 3 || isempty(p), p = struct(); end
if ~isfield(p,'a'), p.a = 2/3; end
if ~isfield(p,'b'), p.b = 4/3; end
if ~isfield(p,'c'), p.c = 1.0; end
if ~isfield(p,'d'), p.d = 1.0; end

x = y(1); z = y(2);
f = [ p.a*x - p.b*x*z;
      p.c*x*z - p.d*z ];
end

