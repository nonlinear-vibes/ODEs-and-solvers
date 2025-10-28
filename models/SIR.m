function f = SIR(~, y)
% SIR  Basic SIR compartmental epidemic model
% Equations:
%        S' = -beta * S * I
%        I' =  beta * S * I - gamma * I
%        R' =  gamma * I

% Parameters
beta  = 0.30;    % infection rate
gamma = 0.10;    % recovery rate

% Unpack states
S = y(1);
I = y(2);
R = y(3);

lambda = beta * S * I;

% System
f  = [ -lambda;
        lambda - gamma * I;
        gamma * I];

end

