function [t, Y] = RKF45(fun, tspan, y0, varargin)
% RK45  Runge–Kutta–Fehlberg 4(5) with adaptive steps.
%
% Sharp tolerances: per-component absolute AND relative tests must both pass.
% The controller adjusts the step size once, and enters an adjustment loop
% only if it needs to be decreased (best by test method).
% Propagated order is selectable:
%   'Order' = 5 (default)  -> accept 5th-order state (standard RKF45)
%   'Order' = 4            -> accept 4th-order state (better for stiff problems)
%
% Usage:
%   [t, y] = RK45(@f, [t0 tf], y0, ...
%                 'AbsTol',1e-6, 'RelTol',1e-3, ...
%                 'InitialStep',1e-3, 'Safety',0.9, ...
%                 'MinStep',1e-12, 'MaxStep',[], ...
%                 'Order',5);
%
% Inputs:
%   fun(t,y)     : RHS function, returns column vector size(y0)
%   tspan        : [t0 tf] with t0 < tf
%   y0           : initial condition (vector)
%
% Options (Name,Value):
%   AbsTol       : scalar or vector (per-component absolute tolerance)
%   RelTol       : scalar relative tolerance
%   InitialStep  : initial step size guess, default (tf - t0)/100
%   Safety       : safety factor in (0,1), default 0.9
%   MinStep      : minimum step size (default 1e-12); if reached and tol fails,
%                  a warning is printed and the step is accepted to proceed.
%   MaxStep      : maximum step size (default tf - t0 if empty)
%   Order        : 4 or 5; which state (y4 or y5) to propagate
%
% Returns:
%   t : [N x 1] times (monotone increasing, ending at tf)
%   Y : [N x ny] states (each row corresponds to t(i))

%% Parse options
p = inputParser;
p.addParameter('AbsTol',      1e-6);
p.addParameter('RelTol',      1e-3);
p.addParameter('InitialStep', []);
p.addParameter('Safety',      0.9);
p.addParameter('MinStep',     1e-12);
p.addParameter('MaxStep',     []);
p.addParameter('Order',       5, @(q) ismember(q,[4,5]));   % which state to keep
p.parse(varargin{:});
opts = p.Results;

%% Basic setup
t0 = tspan(1); tf = tspan(2);
assert(t0 < tf, 'RK45 assumes t0 < tf.');

y  = y0(:);         % column state
ny = numel(y);

% tolerances
AbsTol = opts.AbsTol;
if isscalar(AbsTol), AbsTol = AbsTol*ones(ny,1); else, AbsTol = AbsTol(:); end
assert(numel(AbsTol) == ny, 'AbsTol size must match state dimension.');
RelTol = opts.RelTol;
sfty = opts.Safety;

% initial step and bounds
if isempty(opts.InitialStep), h = (tf - t0)/100; else, h = opts.InitialStep; end
if isempty(opts.MaxStep),     hmax = (tf - t0);  else, hmax = opts.MaxStep;  end
hmin = opts.MinStep;
h    = min(max(h, hmin), hmax);

% storage
cap = 1024; k = 1;
t = zeros(cap,1); t(1) = t0;
Y = zeros(cap,ny); Y(1,:) = y.';

% integration state
tn = t0; yn = y;
pord = 4;              % lower order of the Fehlberg 4(5) pair
minStepWarned = false; % warn only once if stuck at hmin

%% Main loop
while tn < tf

    % clamp step to hit tf; allow final step smaller than hmin
    if tf - tn <= hmin
        h = tf - tn;          % last step may be < hmin
    else
        h = min(h, tf - tn);  % don’t overshoot
        h = max(h, hmin);     % enforce minimum otherwise
        h = min(h, hmax);     % enforce maximum if any
    end

    % first try with current h
    [y4, y5] = rkpair_fehlberg(fun, tn, yn, h);

    % choose which state to advance with
    if opts.Order == 4, y_keep = y4; else, y_keep = y5; end

    % sharp tolerances (per component, both must hold)
    errVec = abs(y5 - y4);
    M      = abs(y5);                              % magnitude per component
    errA   = max(errVec ./ max(AbsTol, eps));      % absolute test
    errR   = max(errVec ./ max(RelTol.*M, eps));   % relative test
    err    = max(errA, errR);                      % ≤1 iff BOTH pass
    if err == 0, err = 1e-16; end                  % avoid division by zero

    % one attempt (order-based, pord=4 -> exponent 1/5)
    facGrow = sfty * err^(-1/(pord+1));
    hTry    = h * facGrow;
    hTry    = min(hTry, tf - tn);  % don't overshoot
    hTry    = max(hTry, hmin);     % enforce minimum otherwise
    hTry    = min(hTry, hmax);     % enforce maximum if any

    if hTry > h * (1 + (1 - sfty)/100) % only if meaningfully increased
        h = hTry;
        [y4, y5] = rkpair_fehlberg(fun, tn, yn, h);
        if opts.Order == 4, y_keep = y4; else, y_keep = y5; end
        errVec = abs(y5 - y4);
        M      = abs(y5);                              % magnitude per component
        errA   = max(errVec ./ max(AbsTol, eps));      % absolute test
        errR   = max(errVec ./ max(RelTol.*M, eps));   % relative test
        err    = max(errA, errR);                      % ≤1 iff BOTH pass
        if err == 0, err = 1e-16; end                  % avoid division by zero
    end

    % shrink until acceptable
    while err > 1
        facGrow  = sfty * err^(-1/(pord+1));
        hNew = h * facGrow;
        hNew = min(hNew, tf - tn);
        hNew = max(hNew, hmin);
        hNew = min(hNew, hmax);


        % no change => we’re stuck (usually at hmin)
        if abs(hNew - h) <= eps
            if ~minStepWarned
                warning('RK45:MinStep', ...
                    'MinStep (h=%.3g) prevents meeting tolerances at t=%.6g; proceeding anyway.', ...
                    hmin, tn);
                minStepWarned = true;
            end
            % accept current trial step (violating tol) to make progress
            break;
        end

        % try with the smaller step
        h = hNew;
        [y4, y5] = rkpair_fehlberg(fun, tn, yn, h);
        
        if opts.Order == 4, y_keep = y4; else, y_keep = y5; end
        errVec = abs(y5 - y4);
        M      = abs(y5);                              % magnitude per component
        errA   = max(errVec ./ max(AbsTol, eps));      % absolute test
        errR   = max(errVec ./ max(RelTol.*M, eps));   % relative test
        err    = max(errA, errR);                      % ≤1 iff BOTH pass
        if err == 0, err = 1e-16; end                  % avoid division by zero
    end

    % accept step
    tn = tn + h;
    yn = y_keep;

    % store
    k = k + 1;
    if k > cap
        cap = round(1.5*cap);
        t = [t; zeros(cap-k+1,1)];
        Y = [Y; zeros(cap-k+1,ny)];
    end
    t(k)   = tn;
    Y(k,:) = yn.';

end

% trim
t = t(1:k);
Y = Y(1:k,:);

end

%% Fehlberg 4(5) embedded pair
function [y4, y5] = rkpair_fehlberg(fun, t, y, h)
% Fehlberg coefficients
c2=1/4; c3=3/8; c4=12/13; c5=1; c6=1/2;
a21 = 1/4;
a31 = 3/32;      a32 = 9/32;
a41 = 1932/2197; a42 = -7200/2197; a43 = 7296/2197;
a51 = 439/216;   a52 = -8;         a53 = 3680/513;   a54 = -845/4104;
a61 = -8/27;     a62 = 2;          a63 = -3544/2565; a64 = 1859/4104; a65 = -11/40;
b4 = [25/216,     0,        1408/2565,   2197/4104,  -1/5,       0];
b5 = [16/135,     0,        6656/12825,  28561/56430, -9/50,     2/55];

k1 = fun(t,           y);
k2 = fun(t + c2*h,    y + h*(a21*k1));
k3 = fun(t + c3*h,    y + h*(a31*k1 + a32*k2));
k4 = fun(t + c4*h,    y + h*(a41*k1 + a42*k2 + a43*k3));
k5 = fun(t + c5*h,    y + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4));
k6 = fun(t + c6*h,    y + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5));

K  = [k1, k2, k3, k4, k5, k6];
y4 = y + h * (K * b4.');
y5 = y + h * (K * b5.');
end