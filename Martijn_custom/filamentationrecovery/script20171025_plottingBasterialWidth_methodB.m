


%%


allDatasetsWidths = [];
allDatasetsLengths = [];
allDatasetsColors =[];
for dataIdx = 1:numel(datasetsPaths)
    
    % simply brute force this
    allWidths=[]; allLengths=[];
    myColors = linspecer(numel(schnitzcells));
    colorsMixed = myColors(randperm(size(myColors,1)),1);
    colorsMixed = myColors(randperm(size(myColors,1)),2);
    colorsMixed = myColors(randperm(size(myColors,1)),3);
    allColors=[];
    
    for schnitzIdx = 1:numel(schnitzcells)

        %%
        lenghts = schnitzcells(schnitzIdx).length_skeleton;
        areas   = schnitzcells(schnitzIdx).area;
        widths  = areas./lenghts;

        allWidths  = [allWidths, widths];
        allLengths = [allLengths, lenghts];
        allColors  = [allColors; repmat(colorsMixed(schnitzIdx,:), [numel(lenghts),1])];
    end

    allDatasetsWidths = [allDatasetsWidths, allWidths];
    allDatasetsLengths = [allDatasetsLengths, allLengths];
    allDatasetsColors = [allDatasetsColors; allColors];
    
end
%%



figure; clf; hold on;
scatter(allDatasetsLengths,allDatasetsWidths,5^2,allDatasetsColors); 
ylim([0,1]);

[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins] = ...
    binnedaveraging({allDatasetsLengths},{allDatasetsWidths},[0:3:35]);

scatter(binCenters,medianValuesForBins,20^2,'k','FaceColor','k')

%%

figure;hold on;
W=.72;plot([[1:30].*W-(W/2).*W]./[1:30],'LineWidth',2);
ylim([0,1]);

