% MATLAB code for J Duan, MM Malakhov, JJ Pellett, IS Phadke, J Barber, JC Blackwood. "Management efficacy in a metapopulation model of white-nose syndrome"

% Our multi-population code follows an object-oriented structure. This is the
% class file for a grid of hibernacula/populations, which provides a "blueprint"
% for creating a grid of Hibernaculum objects and accounts for dispersal between them.

% Instructions: DO NOT RUN THIS FILE -- this file is used by all other multi-population
% program files.
% This file does not require any manual changes, unless you want to change
% the way that dispersal is handled.

% Dependencies: Hibernaculum.m must be in the same directory

% version 08/05/2019
% Copyright (c) 2019 Mykhaylo M. Malakhov

classdef HibernaculumGrid < handle

% HibernaculumGrid is a subclass of handle, which means that all HibernaculumGrid
% objects are handles; the upshot is that copying a HibernaculumGrid object
% will only copy the handle, so all copies will refer to the same object;
% this is necessary so we can create methods that modify properties without
% having to return the object and reassign it to itself at every modification
    
    % declare the properties of HibernaculumGrid;
    % every HibernaculumGrid object contains these two structures
    properties
        value % a Hibernaculum object
        migrationMatrix % a matrix that stores proportions for dispersal
    end
    
    % the methods of HibernaculumGrid: these are functions that act on the
    % HibernaculumGrid object
    methods
        
        % the constructor method that returns a HibernaculumGrid object
        function this = HibernaculumGrid(m,n,beta,phi)
            
            % check to make sure that four arguments were submitted; for
            % some reason MATLAB throws an error if this check is absent
            if nargin == 4
                this = HibernaculumGrid(m,n); % make HibernaculumGrid an m x n grid
                
                % go through the whole grid, creating a Hibernaculum object in each cell
                for i = 1:m
                    for j = 1:n
                        this(i,j).value = Hibernaculum(beta,phi);
                    end
                end
            end
        end
        
        % this method sets up the dispersal matrix, which for convenience
        % I place within the (1,1) HibernaculumGrid object;
        % migrationMatrix is a huge m*n x m*n matrix because we need to
        % account for dispersal from each population to every
        % other population; every dispersal event (where immigration
        % and emigration are two separate events) is represented
        % by a cell in migrationMatrix
        function ResetMigration(this,dispersal)
            
            % we can't access m,n from above, so we calclate them from the
            % size of HibernaculumGrid
            m = size(this,1);
            n = size(this,2);
            
            % preallocate the size of migrationMatrix
            this(1,1).migrationMatrix = zeros(m*n, m*n);
            
            for i = 1:m*n
                for j = 1:m*n
                    
                    % the non-diagonal cells of migrationMatrix represent
                    % the immigration events; we divide the
                    % proportion of bats that leave each hibernaculum
                    % by the total number of hibernacula that they will
                    % move to, and place that value in each non-diagonal cell
                    if i ~= j
                        this(1,1).migrationMatrix(i,j) = dispersal / (m*n - 1);
                    end
                end
            end
            
            % the diagonal cells of migrationMatrix represent the emigration
            % events; we place the proportion of bats that are left
            % in each hibernaculum in the diagonal cells
            for k = 1:m*n
                this(1,1).migrationMatrix(k,k) = 1 - sum(this(1,1).migrationMatrix(:,k));
            end
        end
        
        % this method applies migrationMatrix to the current metapopulation
        function Migration(this)
            
            % we can't access m,n from above, so we calclate them from the
            % size of HibernaculumGrid
            m = size(this,1);
            n = size(this,2);
            
            % this matrix temporarily stores fullPopulationMatrix from
            % every Hibernaculum object together
            tempPopMatrix = zeros(m*n,4);
            
            % copy the current population densities from fullPopulationMatrix
            % from every Hibernaculum object into tempPopMatrix
            i = 1;
            for j = 1:m
                for k = 1:n
                    tempPopMatrix(i,:) = this(j,k).value.fullPopulationMatrix(end,1:4);
                    i = i + 1;
                end
            end
            
            % the columns of tempPopMatrix are the four bat population state
            % variables; multiply each one separately by migrationMatrix
            for i = 1:4
                tempPopMatrix(:,i) = this(1,1).migrationMatrix * tempPopMatrix(:,i);
            end
            
            % make an additional row in fullPopulationMatrix for each
            % Hibernaculum and place the "migrated" population values
            % there, then place the unchanged P population in that row as well;
            % this row will be deleted after running Hibernation in
            % the Hibernaculum class
            i = 1;
            for q = 1:m
                for w = 1:n
                    this(q,w).value.fullPopulationMatrix(end+1,1:4) = tempPopMatrix(i,:);
                    this(q,w).value.fullPopulationMatrix(end,5) = this(q,w).value.fullPopulationMatrix(end-1,5);
                    i = i + 1;
                end
            end
        end
        
        % this method simulates a full year of metapopulation dynamics
        function FullYear(this,yearDays)
            
            % we can't access m,n from above, so we calclate them from the
            % size of HibernaculumGrid
            m = size(this,1);
            n = size(this,2);
            
            % call the PreMigration method on each Hibernaculum object in
            % HibernaculumGrid
            for i = 1:m
                for j = 1:n
                    this(i,j).value.PreMigration(yearDays)
                end
            end
            
            % call the HibernaculumGrid Migration method; note that this
            % method is only called once because it takes care of dispersal
            % for the whole grid
            this.Migration
            
            % now call the PostMigration method on each Hibernaculum object
            % in HibernaculumGrid
            for i = 1:m
                for j = i:n
                    this(i,j).value.PostMigration(yearDays)
                end
            end
        end
        
    end
    
end
