% add files in the 'models' and 'solvers' folders to the path
if ~exist('RK4','file'), run('setup.m'); end

% integration params
tspan = [0, 40];     
h     = 0.05;   

% initial conditions (2D)
y0 = [0, 8];

% Lotka–Volterra
[t1,  y1 ] = IE(@VDP, tspan, y0, h);
[t2,  y2 ] = RK4(@VDP, tspan, y0, h);
[t45, y45] = ode45(@VDP, tspan, y0);   % built-in solver for reference

plots(t1, y1, t2, y2)