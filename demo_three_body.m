%% Three-body demo with RK4, RK45 (Order=4/5), and ode45
% Pick one of the known periodic orbits:
% 1) Tspan = [0, 6.2];   x0 = [1.2; 0; 0; -1.049357510];
% 2) Tspan = [0, 11.2];  x0 = [0.994; 0; 0; -2.0317326295573368];
% 3) Tspan = [0, 17.1];  x0 = [0.994; 0; 0; -2.0015851063790825];

Tspan = [0, 2*17.1];
x0    = [0.994; 0; 0; -2.0015851063790825];

% Tolerances and initial step
AbsTol = 1e-4;
RelTol = 1e-8;
h0     = 1e-3;

% MATLAB ode45
opts45 = odeset('AbsTol',AbsTol,'RelTol',RelTol);
[t_ode45, y_ode45] = ode45(@threeBodyRHS, Tspan, x0, opts45);

% Fixed-step RK4 (for reference)
[t_rk4, y_rk4] = RK4(@threeBodyRHS, Tspan, x0, h0);

% RK45 (Fehlberg 4(5)), accepting 4th-order state
[t_rkf4, y_rkf4] = RKF45(@threeBodyRHS, Tspan, x0, ...
    'AbsTol',AbsTol, 'RelTol',RelTol, 'InitialStep',h0, 'Safety',0.9, 'Order',4);


% XY trajectory overlay
figure('Name','Three-body trajectory comparison'); hold on; grid on;
plot(y_ode45(:,1), y_ode45(:,3), '--', 'DisplayName','ode45');
plot(y_rk4(:,1),   y_rk4(:,3),   '--', 'DisplayName','RK4 (fixed step)');
plot(y_rkf4(:,1),  y_rkf4(:,3),  'k-', 'LineWidth',1.5, 'DisplayName','RKF45 (Order=4)');

% Draw Earth and Moon
mu = 0.012277471; 
nu = 1 - mu;

rEarth = 0.04;
rMoon  = 0.02;
th     = linspace(0, 2*pi, 200);

hEarth = fill(nu + rEarth*cos(th), 0 + rEarth*sin(th), [0.2 0.4 1], ...
              'EdgeColor','none', 'FaceAlpha',0.5, 'DisplayName','Earth');
hMoon  = fill(-mu + rMoon*cos(th), 0 + rMoon*sin(th), [0.8 0.8 0.9], ...
              'EdgeColor','none', 'FaceAlpha',0.9, 'DisplayName','Moon');

legend('Location','best');
xlabel('x'); ylabel('y');
title('Three-body planar trajectory');
axis equal; axis padded;


