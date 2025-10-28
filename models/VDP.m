function dy = VDP(~, y, params)
% VDP  van der Pol oscillator with inline parameter mu.
% State y = [y1; y2].  y1' = y2, y2' = mu*(1 - y1^2)*y2 - y1    

    if nargin < 3 || isempty(params), params = struct(); end
    if ~isfield(params, 'mu'), params.mu = 3.0; end

    x    = y(1);
    xdot = y(2);

    dy = [ xdot;
           params.mu * (1 - x^2) * xdot - x ];
end
