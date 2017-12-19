
function [pdfsForBins, binCenters, dataPerBins]=binnedpdfs(valuesx,valuesy,edgesForX,edgesForY)
% function [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=binnedaveraging(valuesx,valuesy,bins)
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
binnedValues = cell(1,numel(edgesForX)-1);
for i = 1:numel(valuesx)
    % loop over bins
    for binIdx = 1:(numel(edgesForX)-1)
        
        % select data from x,y line that falls into that bin
        
        % find indices of the data in this bin
        applicableIndices = find(valuesx{i}>edgesForX(binIdx) & valuesx{i}<edgesForX(binIdx+1));
        % store all this data in binnedValues
        binnedValues{binIdx} = [binnedValues{binIdx} valuesy{i}(applicableIndices)];        
    end
end
%% Then calculate averages, standard deviation, etc

% loop over bins again
pdfsForBins=NaN(numel(edgesForX)-1,numel(edgesForY)-1);
nanFlag=0;
dataPerBins={};
for i = 1:numel(binnedValues)
    
    % get current binnedValues
    currentBinnedvalues = binnedValues{i};
    
    % filter out NaN values
    nonNanIdxs = ~isnan(currentBinnedvalues);
    currentBinnedvalues = currentBinnedvalues(nonNanIdxs);
    if any(nonNanIdxs)
        nanFlag = nanFlag+1; 
    end
    
    % Calculate pdf 
    pdfsForBins(i,:) = histcounts(currentBinnedvalues,edgesForY);
    
    dataPerBins{i} = currentBinnedvalues;
    
    %{
    % calculate mean for this bin
    meanValuesForBins(i) = mean(currentBinnedvalues);
    % calculate std for this bin
    stdValuesForBins(i) = std(currentBinnedvalues);
    % calculate standard error for this bin (i.e. /sqrt(n))
    stdErrValuesForBins(i) = stdValuesForBins(i)/sqrt(numel(currentBinnedvalues));
    %}   
    
end

% calculate bin centers
binSizes = (edgesForX(2:end)-edgesForX(1:end-1));
binCenters = edgesForX(2:end)-binSizes./2;

% display warning if nan values were filtered out
if nanFlag>0
    warning('Your data contained NaN values, these were ignored..');
end

end