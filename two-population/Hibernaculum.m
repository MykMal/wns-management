% MATLAB code for J Duan, MM Malakhov, JJ Pellett, IS Phadke, J Barber, JC Blackwood. "Management efficacy in a metapopulation model of white-nose syndrome"

% Our multi-population code follows an object-oriented structure. This is the base class
% file, which provides the "blueprint" for a single hibernaculum/population. It uses the
% ode45 ODE solver to simulate the swarming, hibernation, roosting, and birth phases, each
% contained in separate functions below, for a single hibernaculum. This class is invoked by the
% HibernaculumGrid class, which creates a grid of Hibernaculum objects.

% Instructions: DO NOT RUN THIS FILE -- this file is used by all other multi-population
% program files.
% This file does not require any manual changes, unless you want to change
% the base parameters.

% Dependencies: none

% version 07/31/2019
% Copyright (c) 2019 Mykhaylo M. Malakhov

classdef Hibernaculum < handle

% Hibernaculum is a subclass of handle, which means that all Hibernaculum
% objects are handles. The upshot is that copying a Hibernaculum object
% will only copy the handle, so all copies will refer to the same object.
% Doing so allows us to create methods that modify object properties without
% having to return the object and reassign it to itself at every
% modification.
    
    % declare the properties of Hibernaculum;
    % every Hibernaculum object contains these three structures
    properties
        fullPopulationMatrix % a matrix that stores the full population
        fullTimeVector % a vector that stores the full time series
        params % a struct of parameters
    end
    
    % the methods of Hibernaculum: these are functions that act on the
    % Hibernaculum object
    methods
        
        % the constructor method that returns a Hibernaculum object
        function this = Hibernaculum(beta, phi)
            
            % define base parameters:
            % these are all the standard parameters as in MainWNS, but
            % without controls since those are set in a separate method
            
            % vaccinationStrength doesn't exist in MainWNS but is needed here
            % because the PreMigration and PostMigration methods can't
            % access controlIntensities
            this.params.vaccinationStrength = 0;
            
            this.params.kML = 19150;
            this.params.mu = 1 / (8.5 * 365);
            this.params.gamma = (0.5 * 0.95) / 21;
            this.params.tau = (1 / 83);
            this.params.delta = (1 / 60);
            this.params.beta = beta;
            % change the next two lines if using separate phi values for swarming and hibernation
            this.params.phiS = phi;
            this.params.phiH = phi;
            this.params.a = 0.75 / 92;
            this.params.s = 600;
            this.params.kPD = 10^10;
            this.params.eta = 0.5;
            this.params.omega = 50;
            this.params.pdMortality = 0;
            this.params.lambda = 0;
            this.params.epsilon = 1 / (this.params.s * this.params.delta + 1);
        end
        
        % method that resets fullTimeVector and fullPopulationMatrix
        function Reset(this, IC)
            this.fullPopulationMatrix = [IC(1) IC(2) IC(3) IC(4) IC(5)];
            
            % this vector holds the complete time step iterations, and gets
            % updated at the end of every year;
            % note that fullTimeVector(1) corresponds to the initial
            % condition (IC)
            this.fullTimeVector = 0;
        end
        
        % method that sets the controls:
        % since all of our controls are implemented by modifying parameter
        % values, this method simply overwrites the parameters affected by controls
        function SetControl(this, controlIntensities)
            this.params.vaccinationStrength = controlIntensities(5);
            this.params.tau = (1 / 83) * (1 - controlIntensities(3)) * (1 - controlIntensities(4));
            this.params.delta = (1 / 60) * (1 - controlIntensities(2)) * (1 - controlIntensities(4));
            this.params.kPD = 10^10 * (1 - controlIntensities(1));
            this.params.eta = 0.5 * (1 - controlIntensities(3));
            this.params.pdMortality = controlIntensities(4);
            this.params.epsilon = 1 / (this.params.s * this.params.delta + 1);
        end
        
        % note that we force ode45 to return the solution at daily (fixed) time steps,
        % otherwise the time series of the several hibernacula may have different lengths
        
        % the pre-migration method, which contains only the swarming phase
        function PreMigration(this, yearDays)
            
            % call ode45 to run the swarming phase
            [timeVectorSwarming, populationMatrixSwarming] = ode45(@Swarming,yearDays:yearDays + 61,this.fullPopulationMatrix(end,:),[],this.params);
            
            % ode45 returns the initial condition as the first row of the output matrix (and the corresponding first entry of the time vector);
            % thus, appending the output matrix as-is to the end of the initial conditions (as will be done at the end of the method) will result
            % in the initial conditions appearing twice, so I delete the first entry of the output matrix and time vector;
            % this is also why there must be overlap in time when calling ode45: for example, if swarming ends on day 61 then the first day of
            % hibernation (where nothing will happen and the IC will be returned) must be 61 as well
            populationMatrixSwarming(1,:) = [];
            timeVectorSwarming(1) = [];
            
            % append the results onto the end of fullPopulationMatrix and fullTimeVector
            this.fullPopulationMatrix = [this.fullPopulationMatrix; populationMatrixSwarming];
            this.fullTimeVector = [this.fullTimeVector; timeVectorSwarming];
        end
        
        % the post-migration method, with all phases except for swarming
        function PostMigration(this, yearDays)
            
            % call ode45 to run the hibernation phase with output from swarming
            [timeVectorHibernation, populationMatrixHibernation] = ode45(@Hibernation,yearDays + 61:yearDays + 273,this.fullPopulationMatrix(end,:),[],this.params);
            
            % delete doubled values (see note for swarming)
            populationMatrixHibernation(1,:) = [];
            timeVectorHibernation(1) = [];
            
            % the Migration method in HibernaculumGrid, which is always
            % called before PostMigration, adds an additional row to the
            % end of fullPopulationMatrix with the reclassified state
            % variables; that row is used as the IC for Hibernation, and
            % now I delete it since dispersal should happen instantaneously
            this.fullPopulationMatrix(end,:) = [];
            
            % preallocate the size of the roosting IC vector
            roostingIC = zeros(1,5);
            
            % compute the roosting reclassification;
            % these values will get fed into ode45 as the IC for Roosting1
            roostingIC(1) = populationMatrixHibernation(end,1);
            roostingIC(2) = populationMatrixHibernation(end,2) + this.params.epsilon * populationMatrixHibernation(end,3);
            roostingIC(3) = 0;
            roostingIC(4) = populationMatrixHibernation(end,4);
            roostingIC(5) = populationMatrixHibernation(end,5);
    
            % call ode45 to run the (first) regular roosting phase with the
            % reclassified IC
            [timeVectorRoosting1, populationMatrixRoosting1] = ode45(@Roosting,yearDays + 273:yearDays + 309,roostingIC,[],this.params);
            
            % delete the IC (see note for swarming)
            populationMatrixRoosting1(1,:) = [];
            timeVectorRoosting1(1) = [];
            
            % call ode45 to run the birth phase with output from the first regular roosting phase
            [timeVectorBirth, populationMatrixBirth] = ode45(@Birth,yearDays + 309:yearDays + 330,populationMatrixRoosting1(end,:),[],this.params);
            
            % delete doubled values (see note for swarming)
            populationMatrixBirth(1,:) = [];
            timeVectorBirth(1) = [];
            
            % preallocate the size of the vaccination IC vector
            vaccinationIC = zeros(1,5);
            
            % compute the vaccination reclassification
            % these values will get fed into ode45 as the IC for Roosting2
            vaccinationIC(1) = populationMatrixBirth(end,1) * (1 - this.params.vaccinationStrength);
            vaccinationIC(2) = populationMatrixBirth(end,2);
            vaccinationIC(3) = populationMatrixBirth(end,3);
            vaccinationIC(4) = populationMatrixBirth(end,4) + this.params.vaccinationStrength * populationMatrixBirth(end,1);
            vaccinationIC(5) = populationMatrixBirth(end,5);
    
            % call ode45 to run the (second) regular roosting phase with the reclassified IC
            [timeVectorRoosting2, populationMatrixRoosting2] = ode45(@Roosting,yearDays + 330:yearDays + 365,vaccinationIC,[],this.params);
            
            % delete the IC (see note for swarming)
            populationMatrixRoosting2(1,:) = [];
            timeVectorRoosting2(1) = [];
            
            % concatenate all of the phase population matrices into one matrix
            this.fullPopulationMatrix = [this.fullPopulationMatrix; populationMatrixHibernation; populationMatrixRoosting1; populationMatrixBirth; populationMatrixRoosting2];
            
            % concatenate all of the phase time vectors into one vector
            this.fullTimeVector = [this.fullTimeVector; timeVectorHibernation; timeVectorRoosting1; timeVectorBirth; timeVectorRoosting2];
        end
        
        % method that returns the initial total population: this is used
        % for calculating survival percentages
        function initialPop = InitialPopulation(this)
            initialPop = this.fullPopulationMatrix(1,1) + this.fullPopulationMatrix(1,2) + this.fullPopulationMatrix(1,3) + this.fullPopulationMatrix(1,4);
        end
        
        % method that returns the total population size at the current time in the simulation:
        % this is used for calculating survival percentages
        function finalPop = FinalPopulation(this)
            finalPop = this.fullPopulationMatrix(end,1) + this.fullPopulationMatrix(end,2) + this.fullPopulationMatrix(end,3) + this.fullPopulationMatrix(end,4);
        end
        
        % method that returns the full population matrix and time vector at
        % the current time in the simulation
        function [time, pop] = FullPopulation(this)
            pop = this.fullPopulationMatrix;
            time = this.fullTimeVector;
        end
        
    end
    
end


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
