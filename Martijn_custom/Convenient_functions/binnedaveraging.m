
function [meanValuesForBins, binCenters,stdValuesForBins,stderr]=binnedaveraging(valuesy,valuesx,bins)
% function [meanValuesForBins, binCenters,errorValuesForBins]=binnedaveraging(valuesy,valuesx,bins)
%
% INPUTS
% - valuesy{i} gives a series of y values
% - valuesx{i} gives the corresponding x values
%     numel(valuesy) should equal numle(valuesx)
% - bins defines the bins
%
% OUTPUTS
% - meanValuesForBins:  mean value for all datapoints in that bin
% - binCenters:         centers of the bins you gave as input
% - stdValuesForBins:   standard deviation for each bin
% - stdValuesForBins:   standard deviation / sqrt (N) for each bin, where N
%                       is the number of values in that bin
    
%% Organize data per bin

% loop over different x,y lines in data
binnedValues = cell(1,numel(bins)-1);
for i = 1:numel(valuesx)
    % loop over bins
    for binIdx = 1:(numel(bins)-1)
        
        % select data from x,y line that falls into that bin
        
        % find indices of the data in this bin
        applicableIndices = find(valuesx{i}>bins(binIdx) & valuesx{i}<bins(binIdx+1));        
        % store all this data in binnedValues
        binnedValues{binIdx} = [binnedValues{binIdx} valuesy{i}(applicableIndices)];
    end
end
%% Then calculate averages, standard deviation, etc

% loop over bins again
meanValuesForBins=NaN(1,numel(binnedValues));
stdValuesForBins=NaN(1,numel(binnedValues));
stdErrValuesForBins=NaN(1,numel(binnedValues));
for i = 1:numel(binnedValues)
    % calculate mean for this bin
    meanValuesForBins(i) = mean(binnedValues{i});
    % calculate std for this bin
    stdValuesForBins(i) = std(binnedValues{i});
    % calculate standard error for this bin (i.e. /sqrt(n))
    stdErrValuesForBins(i) = stdValuesForBins(i)/sqrt(numel(binnedValues{i}));
end

% calculate bin centers
binCenters=bins(2:end)-(bins(2)-bins(1))/2;

end