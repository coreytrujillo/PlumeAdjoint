function [dJdE] = plume3ADM(x, t, dJ, D, B)

% [dc, dJ] = plume3TLM(x, t, dE, D, B)
%
% This is the adjoint model of plume3.m 
%
% Calculates the differential concentration, dc, and differential cost 
% function, dJ, given a time vector, v, space vector, x, differential 
% emission over time vector, dE, and constants deposition, D and 
% background, B. The output, dc, is a matrix of differential concentrations
% with constant space in the rows and constant time in the columns. 
% J is the sum of the concentrations overa all space at the final time.

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
    
% Final Condition
% c = plume3(x,t,E,D,B);
% dJ = 0;
dc(:,end) = dc(:,end) + dJ;
dJdE = zeros(length(t),1);

for k = nt:-1:2 % Time Loop
    for i = nx:-1:1 % Space loop
%         dc(i,k) = 0.1;
        if i == 1 % Boundary Condition
            dc(i,k-1) = dc(i,k-1)+ ( 1 - sigma )*(1-D)* dc(i,k);
%             dc(i,k) = 0;
        elseif i == floor(nx/4) % Emission location
            dJdE(k) = dJdE(k) + dt*(1-D)*dc(i,k);
            dc(i,k-1) =  dc(i,k-1) + (1-sigma)*(1-D)*dc(i,k);
            dc(i-1,k-1) = dc(i-1,k-1) + sigma*(1-D)*dc(i,k);
%             dc(i,k) = 0;
        else
            dc(i,k-1) = dc(i,k-1) +  (1 - sigma)*(1-D)*dc(i,k);
            dc(i-1,k-1) = dc(i-1,k-1) + sigma*(1-D)*dc(i,k);
%             dc(i,k) = 0;
        end
%         dc(i,k) = 0;
    end

end
% dc(:,1) = B;
% DE = dE;

end
