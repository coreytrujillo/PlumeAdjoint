% This script compares adjoint, tangent linear, and finite difference plume
% models.

% Format
clc; clear; clf; format compact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 0:0.05:10; % Space Vector
t = 0:0.2:20; % Time vector
E = 0.1*(ones(1,floor(length(t)))); % Base Emissions
D = 0.01; % Deposition
B = 0.01; % Background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the forward model
[c, J] = plume3(x,t,E,D,B);

% Iterate Finite Difference and Tangent Linear models through each timestep
for i = 1:length(t)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finite Difference Model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Perturbed Emissions
    p = 0.01; % Perturbation 
    Ep = E; % Base emission
    Ep(i) = E(i)*(p+1); % Perturb Emissions at this timestep

    % Calculate concentrations and cost function based on perturbed
    % emissions
    [cp, Jp] = plume3(x,t,Ep,D,B);
    
    % Calculate dJ/dE for current time step
    dJ = Jp - J;
    dJdE_FD = dJ./(Ep(i)-E(i));
    
    % Plot point at current timestep
    pl_FD = plot(t(i),dJdE_FD,'bo');
    hold on
    title('Finite Difference Model')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tangent Linear Model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set dE_TL
    dE_TL = zeros(size(E)); % Initialize differential emission vector
    dE_TL(i) = 1; % Set current time step to 1
    
    % Get differential concentrations and differential cost function
    [dJdE_TL] = plume3TLM(x,t,dE_TL,D,B);
    
    % Plot point at current timestep
    pl_TL = plot(t(i),dJdE_TL,'kx');
    title('Tangent Linear Model')
    if i==length(t)
        legend('Finite Difference', 'Tangent Linear')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjoint Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate adjoint plume model
[dJdE_AD] = plume3ADM(x,t,1,D,B);

% Plot differnential cost function over entire timespace
pl_AD = plot(t,dJdE_AD, 'm');
legend([pl_FD, pl_TL, pl_AD], 'Finite Difference','Tangent Linear','Adjoint', 'Location', 'northwest')
xlabel('Time')
ylabel('dJ/dE')
title('Comparison of Models')
