% add files in the 'models' and 'solvers' folders to the path
if ~exist('RK4','file'), run('setup.m'); end

% integration params
tspan = [0, 120];     
h     = 0.05;   

% system params
params = struct('mu', 5.0);
f = @(t,y) VDP(t, y, params);

% initial conditions (2D)
y0 = [0, 1];

% Lotka–Volterra
[t1,  y1 ] = EE(f, tspan, y0, h);
[t2,  y2 ] = IE(f, tspan, y0, h);
[t45, y45] = ode45(f, tspan, y0);   % built-in solver for reference

plots(t1, y1, t2, y2)
