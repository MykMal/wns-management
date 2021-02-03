% MATLAB code for J Duan, MM Malakhov, JJ Pellett, IS Phadke, J Barber, JC Blackwood. "Management efficacy in a metapopulation model of white-nose syndrome"

% This is the base function file for single-population WNS simulations.
% This function uses the ode45 ODE solver to simulate the swarming,
% hibernation, roosting and birth phases, each contained in separate
% functions below, for a single hibernaculum.

% Instructions: DO NOT RUN THIS FILE -- this file is used by all other single-subpopulation programs.
% Typically you only want to modify the number of years, but all other model parameters can also be changed here.

% Dependencies: none

% NOTE: The order of controls in controlIntensities here differs from that in the paper and in two-population program files:
% first column is vaccination intensity
% second column is fungicide
% third column is soil bacteria
% fourth column is UV
% fifth column is microclimate

% version 07/31/2019
% Copyright (c) 2019 Mykhaylo M. Malakhov

function [fullTimeVector,populationOutput] = MainWNS(beta,phi,controlIntensities)

years = 2; % number of years to simulate for

S0 = 14999; % initial density of susceptible class
E0 = 1; % initial density of exposed class
I0 = 0; % initial density of infected class
V0 = 0; % initial density of vaccinated class
P0 = 0; % initial density of environmental P. destructans reservoir

% this matrix holds the total population, and gets updated at the end of every year
fullPopulationMatrix = [S0 E0 I0 V0 P0];

% preallocate the size of the roosting initial condition vector
roostingIC = zeros(1,5);

% preallocate the size of the vaccination initial condition vector
vaccinationIC = zeros(1,5);

% this vector holds the complete time step iterations, and gets updated at the end of every year
% note that fullTimeVector(1) corresponds to the initial condition
fullTimeVector = 0;

% define parameters
params.kML = 19150;
params.mu = 1 / (8.5 * 365);
params.gamma = (0.5 * 0.95) / 21;
params.tau = (1 / 83) * (1 - controlIntensities(3)) * (1 - controlIntensities(4));
params.delta = (1 / 60) * (1 - controlIntensities(5)) * (1 - controlIntensities(4));
params.beta = beta;
% change the next two lines if using separate phi values for swarming and hibernation
params.phiS = phi;
params.phiH = phi;
params.a = 0.75 / 92;
params.s = 600;
params.kPD = 10^10 * (1 - controlIntensities(2));
params.eta = 0.5 * (1 - controlIntensities(3));
params.omega = 50;
params.pdMortality = controlIntensities(4);
params.lambda = 0;

% the value of epsilon does not depend on time, so we compute it here
% note that in the paper we write out this function every time it's used instead of calling it epsilon
epsilon = 1 / (params.s * params.delta + 1);

% go through the phases for each year
for i = 1:years
    
    % calculate how many days have already elapsed in past years
    yearDays = (i - 1) * 365;
    
    % call ode45 to run the swarming phase
    [timeVectorSwarming, populationMatrixSwarming] = ode45(@Swarming,[yearDays,yearDays + 61],fullPopulationMatrix(end,:),[],params);
    
    % ode45 returns the initial condition as the first row of the output matrix (and the corresponding first entry of the time vector);
    % thus, appending the output matrix as-is to the end of the initial conditions (as will be done at the end of the year) will result
    % in the initial conditions appearing twice, so we delete the first entry of the output matrix and time vector;
    % this is also why there must be overlap in time when calling ode45: for example, if swarming ends on day 61 then the first day of
    % hibernation (where nothing will happen and the initial condition will be returned) must be 61 as well
    populationMatrixSwarming(1,:) = [];
    timeVectorSwarming(1) = [];
    
    % call ode45 to run the hibernation phase with output from swarming
    [timeVectorHibernation, populationMatrixHibernation] = ode45(@Hibernation,[yearDays + 61,yearDays + 273],populationMatrixSwarming(end,:),[],params);
    
    % delete doubled values (see note for swarming)
    populationMatrixHibernation(1,:) = [];
    timeVectorHibernation(1) = [];
    
    % compute the roosting reclassification
    % these values will get fed into ode45 as the initial condition for Roosting1
    roostingIC(1) = populationMatrixHibernation(end,1);
    roostingIC(2) = populationMatrixHibernation(end,2) + epsilon * populationMatrixHibernation(end,3);
    roostingIC(3) = 0;
    roostingIC(4) = populationMatrixHibernation(end,4);
    roostingIC(5) = populationMatrixHibernation(end,5);
    
    % call ode45 to run the (first) regular roosting phase with the
    % reclassified initial condition
    [timeVectorRoosting1, populationMatrixRoosting1] = ode45(@Roosting,[yearDays + 273,yearDays + 309],roostingIC,[],params);
    
    % delete the initial condition (see note for swarming)
    populationMatrixRoosting1(1,:) = [];
    timeVectorRoosting1(1) = [];
    
    % call ode45 to run the birth phase with output from the first regular roosting phase
    [timeVectorBirth, populationMatrixBirth] = ode45(@Birth,[yearDays + 309,yearDays + 330],populationMatrixRoosting1(end,:),[],params);
    
    % delete doubled values (see note for swarming)
    populationMatrixBirth(1,:) = [];
    timeVectorBirth(1) = [];
    
    % compute the vaccination reclassification;
    % like the roosting reclassification, this does not happen in time
    vaccinationIC(1) = populationMatrixBirth(end,1) * (1 - controlIntensities(1));
    vaccinationIC(2) = populationMatrixBirth(end,2);
    vaccinationIC(3) = populationMatrixBirth(end,3);
    vaccinationIC(4) = populationMatrixBirth(end,4) + controlIntensities(1) * populationMatrixBirth(end,1);
    vaccinationIC(5) = populationMatrixBirth(end,5);
    
    % call ode45 to run the (second) regular roosting phase with the
    % reclassified initial condition
    [timeVectorRoosting2, populationMatrixRoosting2] = ode45(@Roosting,[yearDays + 330,yearDays + 365],vaccinationIC,[],params);
    
    % delete the initial condition (see note for swarming)
    populationMatrixRoosting2(1,:) = [];
    timeVectorRoosting2(1) = [];
    
    % concatenate all of the phase population matrices into one matrix
    fullPopulationMatrix = [fullPopulationMatrix; populationMatrixSwarming; populationMatrixHibernation; populationMatrixRoosting1; populationMatrixBirth; populationMatrixRoosting2];
    % clear the phase population matrices to free up memory
    clearvars populationMatrixSwarming populationMatrixHibernation populationMatrixRoosting1 populationMatrixBirth populationMatrixRoosting2
    
    % concatenate all of the phase time vectors into one vector
    fullTimeVector = [fullTimeVector; timeVectorSwarming; timeVectorHibernation; timeVectorRoosting1; timeVectorBirth; timeVectorRoosting2];
    % clear the phase time vectors to free up memory
    clearvars timeVectorSwarming timeVectorHibernation timeVectorRoosting1 timeVectorBirth timeVectorRoosting2
    
end % ends years loop

% return the survival percentage at the end of time
populationOutput = 100 * ((fullPopulationMatrix(end,1) + fullPopulationMatrix(end,2) + fullPopulationMatrix(end,3) + fullPopulationMatrix(end,4)) / (S0 + E0 + I0 + V0));


% define all of the phase functions that get called by ode45

function swarmingF = Swarming(~,initialSwarmingF,paramsF)
    swarmingF = zeros(5,1); % initialize an empty vector for the five state variables
    
    % here are the swarming phase model equations
    swarmingF(1) = - (paramsF.phiS * initialSwarmingF(5) + paramsF.mu) * initialSwarmingF(1) + paramsF.lambda * initialSwarmingF(4);
    swarmingF(2) = paramsF.phiS * initialSwarmingF(5) * initialSwarmingF(1) - paramsF.mu * initialSwarmingF(2);
    swarmingF(4) = - (paramsF.lambda + paramsF.mu) * initialSwarmingF(4); 
    swarmingF(5) = paramsF.eta * initialSwarmingF(5) * (1 - initialSwarmingF(5) / paramsF.kPD) - paramsF.pdMortality * initialSwarmingF(5);
end

function hibernationF = Hibernation(~,initialHibernationF,paramsF)
    hibernationF = zeros(5,1); % initialize an empty vector for the five state variables
    
    % here are the hibernation phase model equations
    hibernationF(1) = - (paramsF.beta * initialHibernationF(3) + paramsF.phiH * initialHibernationF(5) + paramsF.mu) * initialHibernationF(1) + paramsF.lambda * initialHibernationF(4);
    hibernationF(2) = (paramsF.beta * initialHibernationF(3) + paramsF.phiH * initialHibernationF(5)) * initialHibernationF(1) - (paramsF.tau + paramsF.mu) * initialHibernationF(2);
    hibernationF(3) = paramsF.tau * initialHibernationF(2) - (paramsF.delta + paramsF.mu) * initialHibernationF(3);
    hibernationF(4) = - (paramsF.lambda + paramsF.mu) * initialHibernationF(4);
    hibernationF(5) = (paramsF.omega * initialHibernationF(3) + paramsF.eta * initialHibernationF(5)) * (1 - initialHibernationF(5) / paramsF.kPD) - paramsF.pdMortality * initialHibernationF(5);
end

function roostingF = Roosting(~,initialRoostingF,paramsF)
    roostingF = zeros(5,1); % initialize an empty vector for the five state variables
    
    % here are the regular roosting phase model equations
    roostingF(1) = paramsF.a * initialRoostingF(2) + paramsF.lambda * initialRoostingF(4) - paramsF.mu * initialRoostingF(1);
    roostingF(2) = - (paramsF.a + paramsF.mu) * initialRoostingF(2);
    roostingF(4) = - (paramsF.lambda + paramsF.mu) * initialRoostingF(4);
    roostingF(5) = paramsF.eta * initialRoostingF(5) * (1 - initialRoostingF(5) / paramsF.kPD) - paramsF.pdMortality * initialRoostingF(5);
end

function birthF = Birth(~,initialBirthF,paramsF)
    birthF = zeros(5,1); % initialize an empty vector for the five state variables
    
    % calculate N out here since it is used twice
    NF = initialBirthF(1) + initialBirthF(2) + initialBirthF(3) + initialBirthF(4);
    
    % here are the birth phase model equations
    birthF(1) = paramsF.gamma * NF * (1 - NF / paramsF.kML) + paramsF.a * initialBirthF(2) + paramsF.lambda * initialBirthF(4) - paramsF.mu * initialBirthF(1);
    birthF(2) = - (paramsF.a + paramsF.mu) * initialBirthF(2);
    birthF(4) = - (paramsF.lambda + paramsF.mu) * initialBirthF(4);
    birthF(5) = paramsF.eta * initialBirthF(5) * (1 - initialBirthF(5) / paramsF.kPD) - paramsF.pdMortality * initialBirthF(5);
end

end % ends the main function
