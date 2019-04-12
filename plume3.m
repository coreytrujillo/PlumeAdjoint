function [c, J] = plume3(x, t, E, D, B)

% [c, J] = plume3(x, t, E, D, B)
%
% Calculates the concentration, c, and cost function, J, given a time 
% vector, v, space vector, x, Emission over time vector, E, and constants 
% deposition, D and background, B. The output, c, is a matrix of
% concentrations with constant space in the rows and constant time in the 
% columns. J is the sum of the concentrations overa all space at the final
% time.
% 
% Similar to plume2 but reformatted for ease of use in building TLM and ADM

% Preallocate concentration matrix
c = zeros(length(x), length(t));

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
c(:,1) = B;

for k = 2:nt % Time Loop
    for i = 1:nx % Space loop
        
        if i == 1 % Boundary Condition
            c(i,k) = (( 1 - sigma ) * c(i,k-1) + sigma * B)*(1-D);
            
        elseif i == floor(nx/4) % Emission location
            c(i,k) = ((1 - sigma)*c(i,k-1) + sigma*c(i-1,k-1) + E(k)*dt)*(1-D);
            
        else
            c(i,k) = ((1 - sigma)*c(i,k-1) + sigma*c(i-1,k-1))*(1-D);
        end
    end
end

J = sum(c(:,end)); 

end
