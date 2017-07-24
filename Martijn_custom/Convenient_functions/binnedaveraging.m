
function [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins]=binnedaveraging(valuesx,valuesy,edges)
% function [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=binnedaveraging(valuesx,valuesy,bins)
%
% INPUTS
% - valuesy{i} gives a series of y values
% - valuesx{i} gives the corresponding x values
%     numel(valuesy) should equal numle(valuesx)
% - edges defines the edges of the bins
%
% OUTPUTS
% - meanValuesForBins:  mean value for all datapoints in that bin
% - binCenters:         centers of the bins you gave as input
% - stdValuesForBins:   standard deviation for each bin
% - stdValuesForBins:   standard deviation / sqrt (N) for each bin, where N
%                       is the number of values in that bin
    
%% Organize data per bin

% loop over different x,y lines in data
binnedValues = cell(1,numel(edges)-1);
for i = 1:numel(valuesx)
    % loop over bins
    for binIdx = 1:(numel(edges)-1)
        
        % select data from x,y line that falls into that bin
        
        % find indices of the data in this bin
        applicableIndices = find(valuesx{i}>edges(binIdx) & valuesx{i}<edges(binIdx+1));
        % store all this data in binnedValues
        binnedValues{binIdx} = [binnedValues{binIdx} valuesy{i}(applicableIndices)];        
    end
end
%% Then calculate averages, standard deviation, etc

% loop over bins again
meanValuesForBins=NaN(1,numel(binnedValues));
medianValuesForBins=NaN(1,numel(binnedValues));
stdValuesForBins=NaN(1,numel(binnedValues));
stdErrValuesForBins=NaN(1,numel(binnedValues));
counts=NaN(1,numel(binnedValues));
nanFlag=0;
for i = 1:numel(binnedValues)
    
    % get current binnedValues
    currentBinnedvalues = binnedValues{i};
    
    % filter out NaN values
    nonNanIdxs = ~isnan(currentBinnedvalues);
    currentBinnedvalues = currentBinnedvalues(nonNanIdxs);
    if any(nonNanIdxs)
        nanFlag = nanFlag+1; 
    end
    
    % calculate mean for this bin
    meanValuesForBins(i) = mean(currentBinnedvalues);
    % calculate mean for this bin
    medianValuesForBins(i) = median(currentBinnedvalues);
    % calculate std for this bin
    stdValuesForBins(i) = std(currentBinnedvalues);
    % calculate standard error for this bin (i.e. /sqrt(n))
    stdErrValuesForBins(i) = stdValuesForBins(i)/sqrt(numel(currentBinnedvalues));
    % track counts
    counts(i) = numel(currentBinnedvalues);
end

% calculate bin centers
binSizes = (edges(2:end)-edges(1:end-1));
binCenters = edges(2:end)-binSizes./2;

% display warning if nan values were filtered out
if nanFlag>0
    warning('Your data contained NaN values, these were ignored..');
end

end