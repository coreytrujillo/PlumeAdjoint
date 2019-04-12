function [dJ] = plume3TLM(x, t, dE, D, B)

% [dc, dJ] = plume3TLM(x, t, dE, D, B)
%
% This is the tangent linear model of plume3.m
%
% Calculates the differential concentration, dc, and differential cost 
% function, dJ, given a time vector, v, space vector, x, differential 
% emission over time vector, dE, and constants deposition, D and 
% background, B. The output, dJ is the sum of the differential
% concentrations over all space at the final time.


% Preallocate concentration matrix
dc = zeros(length(x), length(t));

% Constant Velocity
v  = 0.2;

% Differential space and time
DT = diff(t);
dt = DT(1);
DX = diff(x);
dx = DX(1);
nt = length(t);
nx = length(x);

% Lumped constant: Courant limit must be >=1
% See Seinfeld and Pandis Eq 25.131
sigma = v * dt/dx;

if sigma > 1.01
    disp('Sigma is too big!')
elseif sigma < 0
    disp('Sigma is Negative!')
else

end
    
% Initial Condition
% dc(:,1) = B;

for k = 2:nt % Time Loop
    for i = 1:nx % Space loop
        
        if i == 1 % Boundary Condition
            dc(i,k) = ( 1 - sigma )*(1-D)* dc(i,k-1);
            
        elseif i == floor(nx/4) % Emission location
            dc(i,k) = ((1 - sigma)*dc(i,k-1) + sigma*dc(i-1,k-1) + dE(k)*dt)*(1-D);
            
        else
            dc(i,k) = ((1 - sigma)*dc(i,k-1) + sigma*dc(i-1,k-1))*(1-D);
        end
    end
end

dJ = sum(dc(:,end));

end
