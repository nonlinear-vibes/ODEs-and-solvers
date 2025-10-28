function f = Rayleigh(~, y, p)
% Rayleigh oscillator (nonlinear damping on y').
% y1' = y2
% y2' = mu*(1 - y2^2)*y2 - y1
% (Rayleigh form: y'' + y = mu*(1 - (y')^2)*y')

if nargin < 3 || isempty(p), p = struct(); end
if ~isfield(p,'mu'), p.mu = 3.0; end

y1 = y(1); y2 = y(2);
f = [ y2;
      p.mu*(1 - y2^2)*y2 - y1 ];
end
