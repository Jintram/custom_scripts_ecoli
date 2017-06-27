



OUTPUTFOLDER = 'U:\PROJECTS\A_CRPcAMP\figures_matlab\pulsing\';

%{
CUSTOMXLIM=[0,500];
CUSTOMYLIM=[0,1.2];
%}

%%
load('H:\EXPERIMENTAL_DATA_2017\2017-03-22_asc1004_cAMP_pulsing\pos2smallcrop2\data\pos2smallcrop2-Schnitz.mat');

% for comparison I could also load the older data from
%load('G:\EXPERIMENTAL_DATA_2016\2016-12-08_asc990_lac\\pos2crop\data\pos2crop-Schnitz.mat')
% ^ Data from previous experiment with wild type cells on gel pad.

myThreeColors = linspecer(3);

disp('Loading done');

%%

FIELDX='time';
FIELDY='muP9_fitNew_all'

%schnitzcells.
FIELDX='time_atC';
FIELDY='C6_mean'
FIELDY='C2_mean'

FIELDX='time_atY';
FIELDY='Y6_mean'
%FIELDY='Y2_mean'
FIELDY='muP9_fitNew_cycCor'
%size([schnitzcells.muP9_fitNew_cyccor])

TITLES  = {'CRP signal','Const. signal','Growth rate','CRP production rate','Constitutive prod. rate'};
XFIELDS = {'time_atY','time_atC','time_atY','time_atdY','time_atdC'};
YFIELDS = {'Y6_mean','C6_mean','muP9_fitNew_cycCor','dY5_sum','dC5_sum'};

%% plot signal evolving over time for growth, fluor label concentrations and production

% Open figures already here to find them back easily later
figure(1); figure(2); figure(3);
figure(4); figure(5);

outputSimple=[]; theYLim = [];
for fieldIndex = 1:numel(XFIELDS)
    
    FIELDX = XFIELDS{fieldIndex};
    FIELDY = YFIELDS{fieldIndex};
    
    %%
    allX  = [schnitzcells.(FIELDX)]/60;
    allY = [schnitzcells.(FIELDY)];

    %% average line

    [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins] = ...
        binnedaveraging({allX},{allY},0:.1:20);

    %%
    h1=figure(fieldIndex); clf; hold on;
    title(TITLES(fieldIndex));
    
    plot(allX,allY,'.');

    noNanIdx = ~isnan(meanValuesForBins);
    plot(binCenters(noNanIdx),meanValuesForBins(noNanIdx),'-k','LineWidth',3);

    xlim([0,max(allX)]);
    maxAverageValue = max(meanValuesForBins(noNanIdx));
    theYLim{fieldIndex} = [-maxAverageValue/2 maxAverageValue*2];
    ylim(theYLim{fieldIndex});
    %ylim([min(allY),max(allY)*1.3]); 
    %if ~isempty(strfind(FIELDY,'mu'))
    %    ylim([-1,3]); % growht rate ylim
    %end

    xlabel('Time (hrs)');

    MW_makeplotlookbetter(20);

    if exist('OUTPUTFOLDER','var')
        saveas(h1,[OUTPUTFOLDER 'fig' num2str(fieldIndex+1) '.tif'])
    end

    %%
    %{
    %crosscorr(growthRates,CRPsignal)

    % R(growthRates , CRPsignal)
    mySettings=struct;
    [S, NDtau,ccNormalizationFactor]=general_get_crosscorrelation_no_weighing(growthRates,CRPsignal,0,mySettings)
    %}
    
    outputSimple.(FIELDY).x = binCenters(noNanIdx);
    outputSimple.(FIELDY).y = meanValuesForBins(noNanIdx);
    outputSimple.(FIELDY).allx = allX;
    outputSimple.(FIELDY).ally = allY;    
        
end

%% Now produce a block signal of the pulsing

for dataIdx=1:5
    h1=figure(dataIdx);

    EXPERIMENTALSTARTTIME=min([schnitzcells.timestamp])*24;

    cAMPSwitchTimes=...
    {...
    {[2017,03,23,09,29,44],    [1]},... 2100uM
    {[2017,03,23,10,29,46],    [2]},... 43uM
    {[2017,03,23,11,29,49],    [1]},... 2100uM
    {[2017,03,23,12,29,51],    [2]},... 43uM
    {[2017,03,23,13,29,54],    [1]},... 2100uM
    {[2017,03,23,14,29,56],    [2]},... etc
    {[2017,03,23,19,29,59],    [1]},...
    {[2017,03,24,00,30,01],    [2]},...
    {[2017,03,24,05,30,04],    [1]},...
    {[2017,03,24,05,39,40],    [2]},...
    {[2017,03,24,07,58,04],    [1]}...
    };

    switchTimesHrsRaw = arrayfun(@(x) datenum(cAMPSwitchTimes{x}{1}),1:numel(cAMPSwitchTimes))*24;
    switchTimesHrs=switchTimesHrsRaw-EXPERIMENTALSTARTTIME;

    switchTimesHrsCorrected = switchTimesHrs+58/60

    % 

    lineY = theYLim{dataIdx};
    
    %plot(switchTimesHrsCorrected,ones(1,numel(switchTimesHrsCorrected))*100,'ko','MarkerSize',10,'MarkerFaceColor','k');
    lines=[]; myThreeColors = linspecer(3);
    for idx=1:numel(switchTimesHrsCorrected)

        switchColor = myThreeColors(cAMPSwitchTimes{idx}{2}+1,:);

        lines(end+1)=plot([switchTimesHrsCorrected(idx) switchTimesHrsCorrected(idx)],lineY,':o','MarkerSize',10,...
            'MarkerFaceColor',switchColor,'Color',switchColor,'MarkerEdgeColor',switchColor);

    end
    legend(lines(1:2),{'2100uM cAMP','43uM'});
    
    if exist('OUTPUTFOLDER','var')
        saveas(h1,[OUTPUTFOLDER 'fig' num2str(dataIdx+1) '.tif'])
    end
end


%% correlation functions

SETS = {[1,3],[2,3],[1,2]}
SETNAMES = {'R(CRP,\mu)','R(\sigma70,\mu)','R(CRP,\sigma70)'};

for setIdx= 1:numel(SETS)

    % which 2 y-signals do we want cross-corr of?
    currentSet = SETS{setIdx};
    
    % get y data for both signals
    y1 = outputSimple.(YFIELDS{currentSet(1)}).y;
    y2 = outputSimple.(YFIELDS{currentSet(2)}).y;
    
    % calculate cross-correlation
    noise1 = y1-mean(y1);
    noise2 = y2-mean(y2);
    
    MAXLAGS = numel(noise1);

    [Rtau, tau] = xcorr(...
            noise1, ...
            noise2, ...
            MAXLAGS,'coeff');
    
    % plot
    h1=figure(numel(XFIELDS)+1+setIdx); clf; hold on;        
        
    title(SETNAMES{setIdx});
    
    % axes
    plot([-MAXLAGS,MAXLAGS],[0,0],'-k');
    plot([0,0],[-1,1],'-k');

    % correlation function
    plot(tau, Rtau,'-o','LineWidth',3)

    xlim([-MAXLAGS,MAXLAGS]);
    ylim([-1,1]);

    ylabel('Correlation');
    xlabel('Delay');
    MW_makeplotlookbetter(20);
    
    if exist('OUTPUTFOLDER','var')
        saveas(h1,[OUTPUTFOLDER 'fig' num2str(numel(XFIELDS)+1+setIdx) '.tif'])
    end
    
end

%% scatter plots (+binned average)

SETSscatter = {[1,3],[2,3]};

for setIdx= 1:numel(SETSscatter)

    SET = SETSscatter{setIdx};

    % get data
    allY1=outputSimple.(YFIELDS{SET(1)}).ally;
    allY2=outputSimple.(YFIELDS{SET(2)}).ally;

    % plot
    h1=figure(numel(XFIELDS)+1+numel(SETS)+setIdx); clf; hold on;

    MYCOLOR=[0.2157 0.4941 0.7216];
    scatter(allY1,allY2,3^2,'filled',...
                        'MarkerEdgeColor','None',...
                        'MarkerFaceColor',MYCOLOR,...
                        'MarkerFaceAlpha',.2,...
                        'MarkerEdgeAlpha',.2);

    %{
    plot(allY1,allY2,'.');
    %}

    xlabel(TITLES{SET(1)});
    ylabel(TITLES{SET(2)});

    [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts] = ...
            binnedaveraging({allY1},{allY2},linspace(0,max(allY1),40));

    MW_makeplotlookbetter(20);

    %xlim([]);
    ylim([0,2]);

    toShow = counts>50;
    errorbar(binCenters(toShow),meanValuesForBins(toShow),stdValuesForBins(toShow),'k-','LineWidth',2)

	if exist('OUTPUTFOLDER','var')
        saveas(h1,[OUTPUTFOLDER 'fig' num2str(numel(XFIELDS)+1+numel(SETS)+setIdx) '.tif'])
    end
    
end

referenceFigNumber=numel(XFIELDS)+1+numel(SETS)+numel(SETSscatter);

%% Towbin divides CRP signal by constitutive signal 

YFIELDS = {'Y6_mean'    'C6_mean'    'muP9_fitNew_cycCor'}; % see below for order

allCRP=outputSimple.(YFIELDS{1}).ally;
allCONSTI=outputSimple.(YFIELDS{2}).ally;
allGROWTH=outputSimple.(YFIELDS{3}).ally;

normalizedSignals = allCRP./allCONSTI;
desiredBins = linspace(0,max(normalizedSignals),40)

[meanValuesForBinsCRP, binCentersYCRP,stdValuesForBinsCRP,stdErrValuesForBinsCRP, countsCRP] = ...
        binnedaveraging({normalizedSignals},{allGROWTH},desiredBins);
toShow = countsCRP>50;
binDistance=binCentersYCRP(2)-binCentersYCRP(1);

%
h1=figure(101); clf; hold on;
scatter(normalizedSignals,allGROWTH,3^2,'filled',...
                        'MarkerEdgeColor','None',...
                        'MarkerFaceColor',MYCOLOR,...
                        'MarkerFaceAlpha',.2,...
                        'MarkerEdgeAlpha',.2);
errorbar(binCentersYCRP(toShow),meanValuesForBinsCRP(toShow),stdValuesForBinsCRP(toShow),'k-','LineWidth',2)

xlim([min(binCentersYCRP(toShow))-binDistance/2,max(binCentersYCRP(toShow))+binDistance/2]);
ylim([0,2]);

MW_makeplotlookbetter(20);

xlabel('CRP normalized by constitutive');
ylabel('Growth rate (dbl/hr)');

if exist('OUTPUTFOLDER','var')
    saveas(h1,[OUTPUTFOLDER 'fig' num2str(101) '.tif'])
end



%% How does the system evolve over time?

timeValues = outputSimple.('Y6_mean').x;
CRPvalues  = outputSimple.('Y6_mean').y;
growthvalues  = outputSimple.('muP9_fitNew_cycCor').y;

h1=figure(referenceFigNumber+1); clf; hold on; 
% saving the figure
SIZE=[17,6.8]; OFFSET = [2,2]; set(h1,'Units','centimeters','Position',[OFFSET SIZE]*2)
set(h1,'RendererMode','manual','Renderer','Painters');

title('Time evolution');

% identify points where a switch has happened just before
% ===
% First get the switch times that are displayed in the plot
displayIdxs=find(switchTimesHrsCorrected>min(timeValues) & switchTimesHrsCorrected<max(timeValues));
switchTimesHrsCorrectedDisplayed = switchTimesHrsCorrected(displayIdxs);
valuesOfSwitchDisplayed = arrayfun(@(idx) cAMPSwitchTimes{idx}{2}, displayIdxs);
% Then identify the first points next to that switchtime
switchIndices=arrayfun(@(x) find(timeValues>switchTimesHrsCorrectedDisplayed(x),1)-1, 1:numel(switchTimesHrsCorrectedDisplayed))
switchIndices=arrayfun(@(x) find(timeValues>switchTimesHrsCorrectedDisplayed(x),1), 1:numel(switchTimesHrsCorrectedDisplayed))

lastIdx = numel(timeValues)-1;
timeColors = makeColorMap([0 0 1], [0 0 0], lastIdx);
timeColors = colormap(parula(lastIdx));
%arrowPlot(CRPvalues,growthvalues,'number',10,'color',[0 0 0]); hold on;
switchcount=0; switchH = [];
for timeIdx = 1:lastIdx    
    plot([CRPvalues(timeIdx) CRPvalues(timeIdx+1)],[growthvalues(timeIdx) growthvalues(timeIdx+1)],'x-','Color',timeColors(timeIdx,:),'LineWidth',2)
    
    if ismember(timeIdx, switchIndices)
        switchcount=switchcount+1;
        switchColor = myThreeColors(valuesOfSwitchDisplayed(switchcount)+1,:);
        switchH(end+1)=plot(CRPvalues(timeIdx),growthvalues(timeIdx),'s','MarkerSize',15,'LineWidth',3,'Color',switchColor); %'MarkerFaceColor','k',
    end
    
    arrowu=CRPvalues(timeIdx+1) -CRPvalues(timeIdx);
    arrowv=growthvalues(timeIdx+1) -growthvalues(timeIdx);
    %quiver(CRPvalues(timeIdx),growthvalues(timeIdx),arrowu,arrowv);
    %arrow([CRPvalues(timeIdx) growthvalues(timeIdx)],[CRPvalues(timeIdx+1) growthvalues(timeIdx+1)],'Length',15,'Width',2);
end

if ~isempty(switchH) % just for special case when looking at other data
    legend(switchH(1:2),{'2100uM cAMP','43uM'});
end

xlabel('cAMP.CRP reporter');
ylabel('Growth rate');

colormap(timeColors);
hC=colorbar();

ylabel(hC, 'Time (hrs)')

inputSettings.rangeIn = [0,1]; % original range of axis
inputSettings.desiredSpacing = 2; %  desired spacing of axis in new metric
inputSettings.rangeOut = [0,max(timeValues)]; % the desired target range
[tickLocationsOldMetric, correspondingLabels] = labelremapping(inputSettings);

set(hC, 'YTickLabel',correspondingLabels, 'XTick',tickLocationsOldMetric)

MW_makeplotlookbetter(20);

if exist('CUSTOMXLIM','var')
    xlim(CUSTOMXLIM);
end
if exist('CUSTOMYLIM','var')
    ylim(CUSTOMYLIM);
end

if exist('OUTPUTFOLDER','var')
    saveas(h1,[OUTPUTFOLDER 'fig' num2str(referenceFigNumber+1) '.tif'])
end

% Find the closest timepoints at which the switch happened

%%

theidx=find(timeValues>switchTimesHrsCorrected(1),1)
switchTimesHrsCorrected(1)
timeValues(theidx)

%% obtain all traces 
TIMEFIELD='time';

% obtain last schnitzes
lastSchnitzes = MW_getschnitzesinlastframe(schnitzcells);

YFIELDS = {'muP9_fitNew_all','C6_mean_all','Y6_mean_all'}
% obtain traces for all of them
[myTraces] = MW_gettracefromschnitzcellsreverse(schnitzcells,lastSchnitzes,TIMEFIELD,YFIELDS);

%% plot traces
IDX=15;

h1=figure(referenceFigNumber+2); clf; hold on; 
divisionPointsToPlot = find(myTraces(IDX).divisionIndices);
plot(myTraces(IDX).time./60,...
    myTraces(IDX).muP9_fitNew_all,'-','LineWidth',2);
plot(myTraces(IDX).time(divisionPointsToPlot)./60,...
     myTraces(IDX).muP9_fitNew_all(divisionPointsToPlot),'s','LineWidth',3);

xlim([0,max(myTraces(IDX).time./60)]); 
ylim([-1,3]);
 
% plot points in time where we switched medium again
lines=[]; myThreeColors = linspecer(3);
for idx=1:numel(switchTimesHrsCorrected)

    switchColor = myThreeColors(cAMPSwitchTimes{idx}{2}+1,:);

    lines(end+1)=plot([switchTimesHrsCorrected(idx) switchTimesHrsCorrected(idx)],lineY,':o','MarkerSize',10,...
        'MarkerFaceColor',switchColor,'Color',switchColor,'MarkerEdgeColor',switchColor);

end
legend(lines(1:2),{'2100uM cAMP','43uM'});

% cosmetics
MW_makeplotlookbetter(20);
xlabel('Time');
ylabel('Growth rate (dbl/hr)');

saveas(h1,[OUTPUTFOLDER 'fig' num2str(referenceFigNumber+2) '.tif'])

%% now plot time evolution for single cells
IDX=13;

selectionIdxs=~isnan(myTraces(IDX).Y6_mean_all);

time = myTraces(IDX).time(selectionIdxs)./60;
y1   = myTraces(IDX).Y6_mean_all(selectionIdxs);
y2   = myTraces(IDX).muP9_fitNew_all(selectionIdxs);

figNr = referenceFigNumber+2;

highlightPointCategories=valuesOfSwitchDisplayed;

h = timeevolutionplot(time,y1,y2,switchIndices,highlightPointCategories,figNr);

schnitzEnd = myTraces(IDX).lineageSchnitzNrs(end);
title(['Trace that ends in schnitz ' num2str(schnitzEnd)])

% saving the figure
SIZE=[8.8,5.8]; OFFSET = [2,2]; set(h,'Units','centimeters','Position',[OFFSET SIZE]*2);
set(h,'RendererMode','manual','Renderer','Painters');

saveas(h1,[OUTPUTFOLDER 'fig' num2str(figNr) '.tif'])

%% Calculate noise parameters

disp('Some stats:');
y=[schnitzcells(:).muP5_fitNew_all];
selection=~isnan(y);
noiseGrowth = std(y(selection))/mean(y(selection))
meanGrowth     = mean(y(selection))
y=[schnitzcells(:).Y5_mean];
noiseCRP    = std(y)/mean(y)
meanCRP     = mean(y)
y=[schnitzcells(:).C5_mean];
noiseS70    = std(y)/mean(y)
meanS70     = mean(y)

%% Some additional plots 




