
function [meanValuesForBins, binCenters,errorValuesForBins]=binnedaveraging(valuesy,valuesx,bins)
% valuesy{i} gives a series of y values
% valuesx{i} gives the corresponding x values
% bins defines the bins
    
%%

binnedValues = cell(1,numel(bins)-1);
for i = 1:numel(valuesx)
    for binIdx = 1:(numel(bins)-1)
        applicableIndices = find(valuesx{i}>bins(binIdx) & valuesx{i}<bins(binIdx+1));
        
        binnedValues{binIdx} = [binnedValues{binIdx} valuesy{i}(applicableIndices)];
    end
end
%%

meanValuesForBins=NaN(1,numel(binnedValues));
errorValuesForBins=NaN(1,numel(binnedValues));
for i = 1:numel(binnedValues)
    meanValuesForBins(i) = mean(binnedValues{i});
    errorValuesForBins(i) = std(binnedValues{i});
end

binCenters=bins(2:end)-(bins(2)-bins(1))/2;

end