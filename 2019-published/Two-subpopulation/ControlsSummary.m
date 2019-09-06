% MATLAB code for J Duan, MM Malakhov, JJ Pellett, IS Phadke, J Barber, JC Blackwood. "Management Efficacy in a Spatially Dynamic Model of White-nose Syndrome"

% This program creates two color grids with control strategy in subpopulation A
% on the y-axis and control strategy in subpopulation B on the x-axis. Each
% cell is colored to represent the percentage of bats alive at the end of the
% simulation. First, the program creates a HibernaculumGrid object, which in turn
% creates a 1x2 grid of Hibernaculum objects. Then, the program runs simulations
% on the grid twice (for low and high dispersal) for each combination of control
% strategies and displays the results in two subplots.

% Instructions: Set the disease transmission case, low and high dispersal proportions, intervention intensity,
% initial conditions, and number of simulation years below. No changes need to be made to other files.

% Dependencies: Hibernaculum.m, HibernaculumGrid.m, and ContourLine.mat 
% (generated by ContourPlot.m) must be in the same directory.

% version 08/10/2019
% Copyright (c) 2019 Mykhaylo M. Malakhov

clear all

% create a new figure so other open figures don't get overwritten, and do
% not show it because it will be saved automatically
f = figure('visible','off');

% set the disease case
% 1 denotes primarily environment-to-bat transmission
% 2 denotes equal contributions to disease transmission
% 3 denotes primarily bat-to-bat transmission
diseaseCase = 2;

% set the low proportion of bats that leave each hibernaculum
dispersalLow = 0;

% set the high proportion of bats that leave each hibernaculum
dispersalHigh = 0.07;

% set the intervention intensity
controlIntensity = 0.9;

% set the initial condition for subpopulation A
IC1 = [14999 1 0 0 0];

% set the initial condition for subpopulation B
IC2 = [14999 1 0 0 0];

% set the number of years to simulate the model for
years = 10;

% load the 25% survival matrix to later pull out phi, beta pairs for each of the three disease cases
load('ContourLine.mat');

if diseaseCase == 1
    % disease case 1 - primarily environment-to-bat transmission
    beta = contourLine(1, 168);
    phi = contourLine(2,168);
elseif diseaseCase == 2
    % disease case 2 - equal contributions to disease transmission
    beta = contourLine(1,117);
    phi = contourLine(2,117);
else
    % disease case 3 - primarily bat-to-bat transmission
    beta = contourLine(1,57);
    phi = contourLine(2,57);
end

% define the dimensions of HibernaculumGrid
m = 1;
n = 2;

% create a HibernaculumGrid object called grid
grid = HibernaculumGrid(m,n,beta,phi);

% create matrices to store the percent changes in survival for low and high
% dispersal, respectively
gridMatrix1 = zeros(6,6);
gridMatrix2 = zeros(6,6);

% store the two dispersal values we will test
dispersalPercentages = [dispersalLow dispersalHigh];

% loop through control strategies in subpopulation A
for i = 1:6
    
    % reset the controls vector
    control1 = zeros(1,5);
    
    % i = 1 represents no control; update the controls vector for anything
    % other than that
    if i ~= 1
        control1(i-1) = controlIntensity;
    end
    
    % update the parameter values in subpopulation A given the current control
    grid(1,1).value.SetControl(control1)
    
    % loop through control strategies in subpopulation B
    for j = 1:6
        
        % reset the controls vector
        control2 = zeros(1,5);
        
        % i = 1 represents no control; update the controls vector for anything
        % other than that
        if j ~= 1
            control2(j-1) = controlIntensity;
        end
        
        % update the parameter values in subpopulation B given the current control
        grid(1,2).value.SetControl(control2)
        
        % first iteration finds the final population size given low dispersal,
        % and the second does the same for high dispersal
        for k = 1:2
            
            % reset the initial conditions for a new simulation
            grid(1,1).value.Reset(IC1)
            grid(1,2).value.Reset(IC2)
            
            % reset the migration matrix with the current dispersal value
            grid.ResetMigration(dispersalPercentages(k))
            
            % simulate the model for the specified number of years
            for year = 1:years
                
                % as in the single hibernaculum model, this keeps track of how many
                % days have elapsed in past years
                yearDays = (year - 1) * 365;
                
                % call the FullYear method in HibernaculumGrid, which includes all
                % dynamics and migration for a single year
                grid.FullYear(yearDays);
            end
            
            % find the total percentage survived, and place it in the proper
            % entry of the proper gridMatrix
            if k == 1
                gridMatrix1(i,j) = 100 * ((grid(1,1).value.FinalPopulation + grid(1,2).value.FinalPopulation) / (grid(1,1).value.InitialPopulation + grid(1,2).value.InitialPopulation));
            else
                gridMatrix2(i,j) = 100 * ((grid(1,1).value.FinalPopulation + grid(1,2).value.FinalPopulation) / (grid(1,1).value.InitialPopulation + grid(1,2).value.InitialPopulation));
            end
            
        end
    end
end

% create the low dispersal subplot
subplot(1,2,1)

% display gridMatrix1 as an image, scaling it to match the default color map
imagesc([0 5],[0 5],gridMatrix1)
axis square

% make the color map range from 0 to 100
caxis([0 100])

% set the color map to a reversed jet array with 100 colors
colormap(flipud(jet(100)))

% replace the numbers on the axes with letters representing the five controls
xticks([0 1 2 3 4 5])
xticklabels({'N','F','M','B','UV','V'})
yticks([0 1 2 3 4 5])
yticklabels({'N','F','M','B','UV','V'})

% label the axes
xlabel('Control Strategy in Subpopulation B')
ylabel('Control Strategy in Subpopulation A')
title('No Dispersal')

% set subplot font size and axes line width
set(gca, 'FontSize', 30, 'LineWidth', 2);

% create the high dispersal subplot
s2 = subplot(1,2,2);

% save the coordinates of the high dispersal subplot before adding the
% color bar
s2Pos = get(s2,'position');

% display gridMatrix2 as an image, scaling it to match the default color map
imagesc([0 5],[0 5],gridMatrix2)
axis square

% make the color map range from 0 to 100
caxis([0 100])

% set the color map to a reversed jet array with 100 colors
colormap(flipud(jet(100)))

% display a color bar showing the color scale
c = colorbar;

% set the color bar label
c.Label.String = 'Percent Survival';

% adding a color bar changes the size of the plot, so now reset it to the
% pre- color bar size to match the first subplot
set(s2,'position',s2Pos);

% replace the numbers on the axes with letters representing the five controls
xticks([0 1 2 3 4 5])
xticklabels({'N','F','M','B','UV','V'})
yticks([0 1 2 3 4 5])
yticklabels({'N','F','M','B','UV','V'})

% label the axes
xlabel('Control Strategy in Subpopulation B')
ylabel('Control Strategy in Subpopulation A')
title('7% Dispersal')

% set subplot font size and axes line width
set(gca, 'FontSize', 30, 'LineWidth', 2);

% set the size of the saved figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gcf,'PaperPosition', [0, 0, 32, 18]);

% save the figure as a JPG with 450 dpi
print(gcf,'SurvivalGridPlot','-djpeg','-r450')
