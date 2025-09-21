# ODEs-and-solvers

Some classic differential equations and basic numerical solvers

## Contents:
Solvers:

- `EE.m` — Explicit (forward) Euler

- `IEN.m` — Implicit Euler (Newton with numerical Jacobian)

- `TRAP.m` — Implicit trapezoidal rule (Newton with numerical Jacobian)

- `RK4.m` — Classic 4th-order Runge–Kutta

- `RKF45.m `— Runge–Kutta–Fehlberg 4(5) with adaptive steps <br />
$~~~~~~~~~~~~~~~~~~~~~~~~~$ Sharp tolerances: per-component absolute and relative tests must both pass. <br />
$~~~~~~~~~~~~~~~~~~~~~~~~~$ Propagate either the 5th-order state (standard) or the 4th-order state (often a bit more robust on <br />
$~~~~~~~~~~~~~~~~~~~~~~~~~$ tricky trajectories).

Utils:
  
- `NewtonIt.m` — Shared Newton iteration with finite-difference Jacobian

Example ODEs:

- `LV.m` — Lotka–Volterra predator–prey model

- `VDP.m` — van der Pol oscillator

- `Rayleigh.m` — Rayleigh oscillator

- `StiffEqn.m` — simple scalar stiff equation

- `CR3BP.m` — Planar circular restricted three-body problem (Earth–Moon parameters)

Demos:

- `demo.m` — Compare solvers on LV/VDP

- `demo_three_body.m` — RK4/RK45/ode45 on the three-body test


![preview](docs/preview.png)


## Quick start

Run `demo`, it will add the models and solvers to the MATLAB path, call two different solvers for the same system and plots their results.
All solvers expect a right-hand side f(t,y) that returns a column vector the same size as y. Each model keeps its parameters at the **top of the file** for easy tweaking.
You can run your own model by creating `models/MySystem.m`:
```matlab
function f = MySystem(t, y)
% y is a column vector
% Example: y1' = y2, y2' = -y1
f = [y(2); -y(1)];
end
```

Then call any solver:
```matlab
[t, y] = SOLVER(@model, tspan, h, y0, varargin);
```



## Environment:

MATLAB R2018a+ recommended (tested on R2023b).

No toolboxes required beyond base MATLAB.
