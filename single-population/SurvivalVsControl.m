% MATLAB code for J Duan, MM Malakhov, JJ Pellett, IS Phadke, J Barber, JC Blackwood. "Management Efficacy in a Spatially Dynamic Model of White-nose Syndrome"

% This program creates five panels, each of which displays a line plot of population survival percentage
% against control intensity for the three disease transmission cases in the single-population setting.

% Instructions: Set the desired number of years in MainWNS.m before running.
% Nothing needs to be modified here.

% Dependencies: This program requires MainWNS.m and ContourLine.mat to be in the same directory.

% version 08/13/2019
% Copyright (c) 2019 Mykhaylo M. Malakhov

clear all

% create a new figure so other open figures don't get overwritten, and do
% not show it because it will be saved automatically
figure('visible','off')

% load the 25% survival matrix to later pull out phi, beta pairs for each of the three disease cases
load('ContourLine.mat')

% the numbers are positions in ContourLine that represent
% 90%, 50%, and 10% of the phi range, respectively
contourIndices = [168 117 57];

% a cell array storing the line styles for each disease transmission case
lineStyles = {':','--','-'};

% the vector of control intensities
controlSpace = linspace(0,1,100);

% step through each of the five controls, according to the old numbering (see note in MainWNS.m)
for control = 1 : 5
	
	% step through the three disease transmission cases
	for i = 1 : 3
        
        count = 1; % keeps track of the current position in controlSpace
        
        % preallocate the vector that stores population survival for each control intensity
		percentSurvivalVector = zeros(1,100);
		
        % set beta and phi
        beta = contourLine(1,contourIndices(i));
        phi = contourLine(2,contourIndices(i));
        
		% step through the range of intervention intensities
		for j = controlSpace
            
            % initialize the current control strategy at the current intervention intensity
			controlIntensities = zeros(1,5);
			controlIntensities(control) = j;
			
			% call MainWNS to run the simulation
			[~, percentSurvival] = MainWNS(beta,phi,controlIntensities);
			
			% set the percentage survived into the proper place in percentSurvivalVector
			percentSurvivalVector(count) = percentSurvival;
			count = count + 1;
		end
		
		% note that the subplots are arranged in the order of controls as listed in the paper
        switch control
			case 1
				subplot(2,3,6)
				title('Vaccination')
			case 2
				subplot(2,3,2)
				title('Fungicide')
			case 3
				subplot(2,3,4)
				title('Soil Bacteria')
			case 4
				subplot(2,3,5)
				title('UV Light')
			case 5
				subplot(2,3,3)
				title('Microclimate')
        end
        
        % plot the results, using a line style corresponding to the current
        % disease transmission case
		plot(controlSpace, percentSurvivalVector, lineStyles{i}, 'LineWidth', 3)
		hold on
	end
	
	% set axes and make labels
	axis([0 1 0 100])
	xlabel('Intervention Intensity')
	ylabel('Percent Survival')
    
    % set subplot font size and axes line width
    set(gca, 'FontSize', 30, 'LineWidth', 2);
end

% make a legend and place it into the space left blank by the absence of a subplot in the 1,1 position
lgd = legend({'Primarily environment-to-bat', 'Equal contributions', 'Primarily bat-to-bat'},'Position',[0.03 0.65 0.4 0.2]);
lgd.Title.String = 'Disease Transmission';

% set legend font size and axes line width
set(lgd, 'FontSize', 30, 'LineWidth', 2);

% set the size of the saved figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gcf,'PaperPosition', [0, 0, 32, 18]);

% save the figure as a JPG with 450 dpi
print(gcf,'SurvivalVsControl','-djpeg','-r450')
