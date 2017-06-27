function [tickLocationsOldMetric, correspondingLabels] = labelremapping(inputSettings) 
% Input 
% 
% Required:
% inputSettings.rangeIn = [101,1]; 
%       original range of axis
% inputSettings.desiredSpacing = .25;
%       desired spacing of axis in new metric
% Choose, EITHER:
% inputSettings.rangeOut = [0,1];
%       the desired target range
% OR
% inputSettings.conversionFactor = 10;
%       the conversion between old and new units
%
% Use the command 
% set(gca, 'XTick',tickLocationsOldMetric,'XTickLabel',correspondingLabels);
% To set it.
% 
% Example:
% inputSettings.rangeIn = [0,1]; % original range of axis
% inputSettings.desiredSpacing = 2; %  desired spacing of axis in new metric
% inputSettings.rangeOut = [0,max(timeValues)]; % the desired target range
% [tickLocationsOldMetric, correspondingLabels] = labelremapping(inputSettings) 

%% 
if ~exist('inputSettings','var')
    inputSettings.rangeIn = [101,1];
    inputSettings.rangeOut = [0,1];
    inputSettings.desiredSpacing = .25;
    inputSettings.conversionFactor = 10;
end

if ~isfield(inputSettings,'desiredSpacing')
    error('Please set inputSettings.desiredSpacing.');
end
if ~isfield(inputSettings,'desiredDecimalsTicks')
    inputSettings.desiredDecimalsTicks=1;
end

if isfield(inputSettings,'rangeOut')
            
    disp('Calculating conversion factor.');
    conversionFactor = (max(inputSettings.rangeOut)-min(inputSettings.rangeOut))/(max(inputSettings.rangeIn)-min(inputSettings.rangeIn));
    
    topValue    = inputSettings.rangeOut(2);
    bottomValue = inputSettings.rangeOut(1);
    distance    = topValue-bottomValue;
    
    tickLocationsOldMetric  = linspace(inputSettings.rangeIn(1),inputSettings.rangeIn(2),distance/(inputSettings.desiredSpacing)+1);
    %tickLocationsNewMetric = linspace(bottomValue,topValue,distance/(inputSettings.desiredSpacing)+1)
    tickLocationsNewMetric = bottomValue:inputSettings.desiredSpacing:topValue;
    
    correspondingLabels = arrayfun(@(x) sprintf(['%0.' num2str(inputSettings.desiredDecimalsTicks) 'f'],tickLocationsNewMetric(x)), 1:numel(tickLocationsNewMetric) ,'UniformOutput',0);
        
    % Only monotic increase allowed, so if not, reverse
    if tickLocationsOldMetric(end)<tickLocationsOldMetric(1)
        tickLocationsOldMetric = tickLocationsOldMetric(end:-1:1);
        correspondingLabels = {correspondingLabels{end:-1:1}};
    end
    
elseif isfield(inputSettings,'conversionFactor')
    
    topValue    = inputSettings.rangeIn(2)*inputSettings.conversionFactor;
    bottomValue = inputSettings.rangeIn(1)*inputSettings.conversionFactor;
    
    correspondingLabels = 0;
    
end

% tickLocationsOldMetric, correspdongingLabels

end