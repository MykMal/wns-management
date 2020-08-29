% MATLAB code for J Duan, MM Malakhov, JJ Pellett, IS Phadke, J Barber, JC Blackwood. "Management Efficacy in a Metapopulation Model of White-nose Syndrome"

% Run this program file to make a heat map of population survival with varying phi, beta for the single-population case.
% This program file also generates ContourLine.mat, which contains the 25% survival contour values.

% Instructions: Set the desired number of years in MainWNS.m before running.
% Nothing needs to be modified here.

% Dependencies: This program requires MainWNS.m to be in the same directory.

% version 07/27/2019
% Copyright (c) 2019 Mykhaylo M. Malakhov

clear all

% create a new figure so other open figures don't get overwritten, and do
% not show it because it will be saved automatically
figure('visible','off')

n = 100; % how many values of beta and phi to iterate through; 100 is what we used

controlIntensities = [0 0 0 0 0]; % keep this at 0 when generating ContourLine.mat

outputXYAxes = zeros(n+1,2); % the zeros will be replaced with beta (x-axis) and phi (y-axis)
contourMatrix = zeros(n+1,n+1); % the output matrix that will be fed into contour

betaMAX = 2 * 10^(-5); % max value of beta to iterate through
phiMAX = 0.8 * 10^(-12); % max value of phi to iterate through

% the two iterators (iteratorY and iteratorX) are for keeping track of what for loop iteration we're currently in
iteratorY = 1; % phi loop iterations

% loop through all phi values
for phi = 0:(phiMAX/n):phiMAX
    
    iteratorX = 1; % beta loop iterations

	% loop through all beta values
    for beta = 0:(betaMAX/n):betaMAX
        
        % MainWNS returns [fullTimeVector,populationOutput] so the ~ simply avoids accepting the time vector
        [~,percentSurvival] = MainWNS(beta,phi,controlIntensities);
        
        % here x and y need to be reversed; it is equivalent to transposing
        % contourMatrix, which is necessary because of how contour works
        contourMatrix(iteratorY,iteratorX) = percentSurvival;
        
        outputXYAxes(iteratorX,1) = beta; % sets first column to be the beta values
        iteratorX = iteratorX + 1;
        
    end % ends beta loop
    
    outputXYAxes(iteratorY,2) = phi; % sets second column to be the phi values
    iteratorY = iteratorY + 1;
    
end % ends phi loop

% plot the contour map of phi vs beta with no lines between regions
contour(outputXYAxes(:,1),outputXYAxes(:,2),contourMatrix,10000,'Fill','on','LineColor','flat');
hold on
% now plot the 25% survival line
[contourLine, ~] = contour(outputXYAxes(:,1),outputXYAxes(:,2),contourMatrix,[25 25],'LineColor','k','LineWidth',4);

% add a color bar and label it
c = colorbar;
caxis([0 100])
c.Label.String = 'Percent Survival';

% delete the first column of contourLine since it has weird values that we are ignoring
contourLine(:,1) = [];

% save contourLine to a matrix file to use in other programs
save('ContourLine.mat','contourLine')

% plot the three disease case points
for i = [168 117 57]
	plot(contourLine(1,i),contourLine(2,i),'h','MarkerSize',24,'MarkerEdgeColor','white','MarkerFaceColor','white')
end

% preserve correct plot aspect ratio
axis image

% label axes
xlabel('\beta')
ylabel('\phi')

% set figure font size and axes line width
set(gca, 'FontSize', 40, 'LineWidth', 2.5);

% set the size of the saved figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gcf,'PaperPosition', [0, 0, 20, 20]);

% save the figure as a JPG with 450 dpi
print(gcf,'25SurvivalContour','-djpeg','-r450')
