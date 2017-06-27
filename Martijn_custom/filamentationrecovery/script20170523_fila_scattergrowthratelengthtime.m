


%%
TRESHOLDTIMEPOINT = 500;

myColorMap = makeColorMap([1 0 0],[.7 .7 .7],[0 0 1],100); % grey-blue colormap

% figure
if exist('FIGURENR','var')
    hScatterLenMuTime=figure(FIGURENR); clf; hold on;
else
    hScatterLenMuTime=figure; clf; hold on;
end

%
% loop over datasets
for dataSetIndex = 1:numel(datasetsPaths)

    %%
    clear schnitzcells
    if ~SPECIALCASE
        schnitzcells = loadandrename(datasetsPaths{dataSetIndex});        
    else
        schnitzcells=S_all_shifted{dataSetIndex};
    end

    
    %% get data
    
    currentSwitchTime = switchTimes(dataSetIndex);
    
    allTimeData = [schnitzcells(:).time]-currentSwitchTime;
    timePoints = unique(allTimeData);
    allLenData = [schnitzcells(:).(LENGTHFIELD)];
    allGrowthRateData = [schnitzcells(:).(GROWTHRATEFIELD)];

    %% determine means (not plotted any more)
    %{
    meanLength = NaN(1,numel(timePoints));
    meanGrowthRate = NaN(1,numel(timePoints));
    for idx = 1:numel(timePoints)

        % brute force get which datapoints belong to this time
        selectedDataIdxs = allTimeData==timePoints(idx);

        % determine mean length of cells in colony
        meanLength(idx) = mean(allLenData(selectedDataIdxs));

        % determine mean growth rate of cells in colony
        meanGrowthRate(idx) = mean(allGrowthRateData(selectedDataIdxs));
    end

    %% plot of how it evolves

    figure; clf; hold on;
    yyaxis left
    plot(timePoints,meanLength,'LineWidth',2);
    MW_makeplotlookbetter(20);
    yyaxis right
    plot(timePoints,meanGrowthRate,'LineWidth',2);
    MW_makeplotlookbetter(20);
    
    xlim([0,max(allTimeData)]);
    %}
    
    %% scatter plot of growth rate vs. length, color-coded by time from switch
    
    selectedDataPoints = allTimeData>0;
    
    % color-coding for time datapoints
    timePointsColored = allTimeData(selectedDataPoints);
    timePointsColored = timePointsColored./TRESHOLDTIMEPOINT;
    timePointsColored(timePointsColored>1) = 1;
    timePointsColored = round(timePointsColored*99)+1;
    
    % plot scatter
    figure(hScatterLenMuTime);     
    scatter(allLenData(selectedDataPoints), allGrowthRateData(selectedDataPoints),2^2,myColorMap(timePointsColored,:),'filled');

    % labels, cosmetics
    ylim([0,2]);
    xlim([0,50]);
    MW_makeplotlookbetter(20);
    xlabel('Cell length (\mum)');
    ylabel('Growth rate (dbl/hr)'); 
   
end
   

% color bar
hTheColorbarWithScatter = figure;
colormap(myColorMap)    
hcb=colorbar;
max(timePointsColored)    
inputSettings.rangeIn = [0,1];  % original range of axis
inputSettings.desiredSpacing = 2; % desired spacing of axis in new metric  
inputSettings.rangeOut = [0,TRESHOLDTIMEPOINT]/60; %         
[tickLocationsOldMetric, correspondingLabels] = labelremapping(inputSettings); % custom function
set(hcb,'YTick',tickLocationsOldMetric,'YTickLabel',correspondingLabels)
hcb.Label.String = 'Time after stress relieve (hrs)';
MW_makeplotlookbetter(20);



%% 