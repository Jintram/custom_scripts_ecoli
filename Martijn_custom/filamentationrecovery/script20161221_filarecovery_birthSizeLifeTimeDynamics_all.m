
NOSAVEPLEASE=1;
YLIMbirthlife=[0,200];

myPlotting.style = 'CBmanuscript';

%% dataset I
WHATDATA = 'tetracycline';
LENGTHFIELD = 'length_fitNew';

RUNSECTIONSFILADIV = 'loadData';
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'lifeTimeVsBirthLength2';
script20160429_filamentRecoveryDivisionRatioss

combinedDynamicsData.(WHATDATA).binCentersDynData = binCentersDynData;
combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData = meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData = stdValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdErrValuesForBins = stdErrValuesForBins;
combinedDynamicsData.(WHATDATA).CVForBinsDynData = stdValuesForBinsDynData./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).VarByMeanForBinsDynData = stdValuesForBinsDynData.^2./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).selectedXdata = selectedXdata; % selected for after switchtime
combinedDynamicsData.(WHATDATA).selectedYdata = selectedYdata;

% figure;plot(combinedDynamicsData.(WHATDATA).binCentersDynData,combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData.*combinedDynamicsData.(WHATDATA).binCentersDynData,'o'); ylim([0,max(combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData.*combinedDynamicsData.(WHATDATA).binCentersDynData)*1.2]);

%% dataset II
WHATDATA = 'temperature';
LENGTHFIELD = 'length_skeleton';

RUNSECTIONSFILADIV = 'loadData';
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'lifeTimeVsBirthLength2';
script20160429_filamentRecoveryDivisionRatioss

combinedDynamicsData.(WHATDATA).binCentersDynData = binCentersDynData;
combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData = meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData = stdValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdErrValuesForBins = stdErrValuesForBins;
combinedDynamicsData.(WHATDATA).CVForBinsDynData = stdValuesForBinsDynData./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).VarByMeanForBinsDynData = stdValuesForBinsDynData.^2./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).selectedXdata = selectedXdata; % selected for after switchtime
combinedDynamicsData.(WHATDATA).selectedYdata = selectedYdata;

%% dataset III
WHATDATA = 'sulA';
LENGTHFIELD = 'length_skeleton';

RUNSECTIONSFILADIV = 'loadData';
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'lifeTimeVsBirthLength2';
script20160429_filamentRecoveryDivisionRatioss

combinedDynamicsData.(WHATDATA).binCentersDynData = binCentersDynData;
combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData = meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData = stdValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdErrValuesForBins = stdErrValuesForBins;
combinedDynamicsData.(WHATDATA).CVForBinsDynData = stdValuesForBinsDynData./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).VarByMeanForBinsDynData = stdValuesForBinsDynData.^2./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).selectedXdata = selectedXdata; % selected for after switchtime
combinedDynamicsData.(WHATDATA).selectedYdata = selectedYdata;

%% dataset IV
WHATDATA = 'deltaMinTET';
LENGTHFIELD = 'length_skeleton';

RUNSECTIONSFILADIV = 'loadData';
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'lifeTimeVsBirthLength2';
script20160429_filamentRecoveryDivisionRatioss

combinedDynamicsData.(WHATDATA).binCentersDynData = binCentersDynData;
combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData = meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData = stdValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdErrValuesForBins = stdErrValuesForBins;
combinedDynamicsData.(WHATDATA).CVForBinsDynData = stdValuesForBinsDynData./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).VarByMeanForBinsDynData = stdValuesForBinsDynData.^2./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).selectedXdata = selectedXdata; % selected for after switchtime
combinedDynamicsData.(WHATDATA).selectedYdata = selectedYdata;

%% Data by Donachie

donachieInterdivData=[...
3.997882242852569, 72.84899899037158;...
4.968603019035189, 34.26358688960575;...
6.126474426851191, 41.085350538058066;...
7.12871530941417, 41.09483119505528;...
7.927799256322489, 23.12941712428279;...
9.00490039153882, 27.30176808096727;...
9.920460981555811, 26.175293654116082;...
11.101233716663794, 27.132409071880616;...
12.12711467901204, 21.088059297200118;...
13.1448693639341, 17.12471373340884;...
14.041961141618854, 20.72779433130586;...
15.065625846487228, 15.25099115959516;...
16.164150804008962, 13.937058287571716;...
16.967913516708116, 18.106823610529688;...
18.239799059321818, 19.06480090620306;...
20.027333842251714, 27.973601910906446;...
19.079760644192177, 13.96463838065452;...
21.087936171784577, 13.037688689699323;...
22.090177054347556, 13.047169346696535;...
23.26356226452264, 15.896106774360348;...
24.0028072594745, 13.254450983772074;...
26.09618557462631, 13.841820778645129];

%% Old version of graph

hOld=figure; clf; hold on; 

myColors = linspecer(5);

DATATOPLOT = {'tetracycline', 'temperature', 'sulA'}%,'deltaMinTET'}
%DATATOPLOT = {'tetracycline'}

legendLines=[];
% Plot large circle + error bar if multiple datapoints
for dataSetindex = 1:numel(DATATOPLOT)
    
    % Parameter for dataset selection
    WHATDATA=DATATOPLOT{dataSetindex}
    % And color
    theColor = myColors(dataSetindex,:);
    
    % Select data that have multiple points
    pointsWithMultipleDataIndices = combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData>0
    
    % Plots those
    xdata   = combinedDynamicsData.(WHATDATA).binCentersDynData(pointsWithMultipleDataIndices);
    ydata   = combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData(pointsWithMultipleDataIndices);
    errdata = combinedDynamicsData.(WHATDATA).stdErrValuesForBins(pointsWithMultipleDataIndices);
    errorbar(xdata,ydata,errdata,...
            '-o','Color',theColor,'MarkerSize',5,'LineWidth',2);%,'LineWidth',1)
    %plot(combinedDynamicsData.(WHATDATA).binCentersDynData,combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData,'.','Color',theColor,'LineWidth',2);
    combinedDynamicsData.(WHATDATA).binCentersDynData(1)
end
% single data on top
for dataSetindex = 1:numel(DATATOPLOT)
    
    WHATDATA=DATATOPLOT{dataSetindex};
    theColor = myColors(dataSetindex,:);    
    
    %legendLines(end+1)=plot(combinedDynamicsData.(WHATDATA).binCentersDynData,combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData,'x','Color',theColor,'LineWidth',3,'MarkerSize',10);    
    
    selectedXdata = combinedDynamicsData.(WHATDATA).selectedXdata;
    selectedYdata = combinedDynamicsData.(WHATDATA).selectedYdata;
    %scatter(selectedXdata,selectedYdata,'.','Color',theColor); % selected for after switchtime

end

ylim([0,100]);
xlabel('Birth size (µm)');
ylabel('Interdivision time (min)');
MW_makeplotlookbetter(20);

legend(legendLines,{'Tetracycline', 'Temperature', 'SulA'});
legend boxoff  

%hFigSIIC = h;

%% Additional statistical analysis
if exist('ALSOPLOTDELTAMINFLAG2')
    ALSOPLOTDELTAMIN=ALSOPLOTDELTAMINFLAG2;
else
    ALSOPLOTDELTAMIN=1;
end
h2=figure(2); clf; hold on;

NRBINSFORSTATS=20;
POINTSTRESHOLD=10; % minimum nr of points to create violin plot

disp(['In violin plot comparison, points with less than ' num2str(POINTSTRESHOLD) ' points are shown separately.']);

minX=min([combinedDynamicsData.('tetracycline').selectedXdata,combinedDynamicsData.('temperature').selectedXdata,combinedDynamicsData.('sulA').selectedXdata]);
maxX=max([combinedDynamicsData.('tetracycline').selectedXdata,combinedDynamicsData.('temperature').selectedXdata,combinedDynamicsData.('sulA').selectedXdata]);
edges=linspace(0,40,NRBINSFORSTATS+1);

% Get groups of data according to bins
%[pdfsForBinsTET, binCentersTET, dataPerBinsTET]=binnedpdfs({combinedDynamicsData.('tetracycline').selectedXdata},{combinedDynamicsData.('tetracycline').selectedYdata},edges,[0,inf]);
%[pdfsForBinsTemp, binCentersTemp, dataPerBinsTemp]=binnedpdfs({combinedDynamicsData.('temperature').selectedXdata},{combinedDynamicsData.('temperature').selectedYdata},edges,[0,inf]);
%[pdfsForBinsSUL, binCentersSUL, dataPerBinsSUL]=binnedpdfs({combinedDynamicsData.('sulA').selectedXdata},{combinedDynamicsData.('sulA').selectedYdata},edges,[0,inf]);
[meanValuesForBinsTET, binCentersTET,stdValuesForBinsTET,stdErrValuesForBinsTET, countsTET,dataPerBinsTET,medianValuesForBinsTET]=binnedaveraging({combinedDynamicsData.('tetracycline').selectedXdata},{combinedDynamicsData.('tetracycline').selectedYdata},edges);
[meanValuesForBinsTemp, binCentersTemp,stdValuesForBinsTemp,stdErrValuesForBinsTemp, countsTemp,dataPerBinsTemp,medianValuesForBinsTemp]=binnedaveraging({combinedDynamicsData.('temperature').selectedXdata},{combinedDynamicsData.('temperature').selectedYdata},edges);
[meanValuesForBinsSUL, binCentersSUL,stdValuesForBinsSUL,stdErrValuesForBinsSUL, countsSUL,dataPerBinsSUL,medianValuesForBinsSUL]=binnedaveraging({combinedDynamicsData.('sulA').selectedXdata},{combinedDynamicsData.('sulA').selectedYdata},edges);
[meanValuesForBinsDMin, binCentersDMin,stdValuesForBinsDMin,stdErrValuesForBinsDMin, countsDMin,dataPerBinsDMin,medianValuesForBinsDMin]=binnedaveraging({combinedDynamicsData.('deltaMinTET').selectedXdata},{combinedDynamicsData.('deltaMinTET').selectedYdata},edges);

% Plotting of raw data and distributions
MYCOLORS=linspecer(5);

subplot(2,1,1); hold on;
l1=plot(combinedDynamicsData.('tetracycline').selectedXdata,combinedDynamicsData.('tetracycline').selectedYdata,'.','Color',MYCOLORS(1,:));
l2=plot(combinedDynamicsData.('temperature').selectedXdata,combinedDynamicsData.('temperature').selectedYdata,'.','Color',MYCOLORS(2,:));
l3=plot(combinedDynamicsData.('sulA').selectedXdata,combinedDynamicsData.('sulA').selectedYdata,'.','Color',MYCOLORS(3,:));
if exist('ALSOPLOTDELTAMIN','var')
    l4=plot(combinedDynamicsData.('deltaMinTET').selectedXdata,combinedDynamicsData.('deltaMinTET').selectedYdata,'.','Color',MYCOLORS(4,:));
end

xlim([0,20]);

if ~exist('ALSOPLOTDELTAMIN','var')
    legend({'Tetracycline','Temperature','SulA'});
else
    legend({'Tetracycline','Temperature','SulA','Min mutant'});
end
        

MW_makeplotlookbetter(16);

% Violin plots
% ===
subplot(2,1,2); hold on;
%figure(3); clf; hold on;
%lala=arrayfun(@(i)dataPerBinsTET{i}',1:numel(dataPerBinsTET),'UniformOutput',0)

% Tetracycline
enoughPointsTET=[arrayfun(@(i) numel(dataPerBinsTET{i})>POINTSTRESHOLD, 1:numel(dataPerBinsTET))];
violin({dataPerBinsTET{enoughPointsTET}},'x',binCentersTET(enoughPointsTET),'facecolor','r','edgecolor',MYCOLORS(1,:),'facealpha',0,'mc',[],'medc',[]);
for i=find(~enoughPointsTET) % plot single points for bins w. little points
    plot(binCentersTET(i)*ones(1,numel(dataPerBinsTET{i})),dataPerBinsTET{i},'.','Color',MYCOLORS(1,:));
end
% Temperature
enoughPointsTemp=[arrayfun(@(i) numel(dataPerBinsTemp{i})>POINTSTRESHOLD, 1:numel(dataPerBinsTemp))];
violin({dataPerBinsTemp{enoughPointsTemp}},'x',binCentersTemp(enoughPointsTemp),'facecolor','g','edgecolor',MYCOLORS(2,:),'facealpha',0,'mc',[],'medc',[]);
for i=find(~enoughPointsTemp) % plot single points for bins w. little points
    plot(binCentersTemp(i)*ones(1,numel(dataPerBinsTemp{i})),dataPerBinsTemp{i},'.','Color',MYCOLORS(2,:));
end
% SulA
enoughPointsSUL=[arrayfun(@(i) numel(dataPerBinsSUL{i})>POINTSTRESHOLD, 1:numel(dataPerBinsSUL))];
violin({dataPerBinsSUL{enoughPointsSUL}},'x',binCentersSUL(enoughPointsSUL),'facecolor','b','edgecolor',MYCOLORS(3,:),'facealpha',0,'mc',[],'medc',[]);
for i=find(~enoughPointsSUL) % plot single points for bins w. little points
    plot(binCentersSUL(i)*ones(1,numel(dataPerBinsSUL{i})),dataPerBinsSUL{i},'.','Color',MYCOLORS(3,:));
end
% DMin
enoughPointsDMin=[arrayfun(@(i) numel(dataPerBinsDMin{i})>POINTSTRESHOLD, 1:numel(dataPerBinsDMin))];
violin({dataPerBinsDMin{enoughPointsDMin}},'x',binCentersDMin(enoughPointsDMin),'facecolor','b','edgecolor',MYCOLORS(4,:),'facealpha',0,'mc',[],'medc',[]);
for i=find(~enoughPointsDMin) % plot single points for bins w. little points
    plot(binCentersDMin(i)*ones(1,numel(dataPerBinsDMin{i})),dataPerBinsDMin{i},'.','Color',MYCOLORS(4,:));
end

xlim([0,20]);
ylim([0,200]);
%violion([combinedDynamicsData.('tetracycline').selectedYdata',...
%   combinedDynamicsData.('tetracycline').selectedYdata']);

MW_makeplotlookbetter(16);

ylabel('Interdivision time (minutes)');
xlabel('Birth size (µm)');

if ~exist('ALSOPLOTDELTAMIN','var')
    legend([l1 l2 l3],{'Tetracycline','Temperature','SulA'});
else
    legend([l1 l2 l3,l4],{'Tetracycline','Temperature','SulA','Min mutant'});
end

% Perform & plot tests
% ===

% Perform t-testing
% Note that bins are equal
% Comparing TET with Temp
for i=find(enoughPointsTET&enoughPointsTemp)
    outcomeNotSameTemp=ttest2(dataPerBinsTET{i},dataPerBinsTemp{i});
    if outcomeNotSameTemp
        plot(binCentersTemp(i),190,'*','Color',MYCOLORS(2,:));
    end
end
% Comparing TET with Sul
for i=find(enoughPointsTET&enoughPointsSUL)
    outcomeNotSameSUL=ttest2(dataPerBinsTET{i},dataPerBinsSUL{i});
    if outcomeNotSameSUL
        plot(binCentersSUL(i),180,'*','Color',MYCOLORS(3,:));
    end
end
if exist('ALSOPLOTDELTAMIN','var')
% Comparing TET with DMin
for i=find(enoughPointsTET&enoughPointsDMin)
    outcomeNotSameDMin=ttest2(dataPerBinsTET{i},dataPerBinsDMin{i});
    if outcomeNotSameDMin
        plot(binCentersDMin(i),170,'*','Color',MYCOLORS(4,:));
    end
end
end

hStatisticsTimeCollapse=h2;

%% Newer version of graph
clear ALSOPLOTDELTAMIN;

h=figure(4); clf; hold on; 

% Actual errorbar figure
moreThanOnePoint=countsTET>1;
errorbar(binCentersTET(moreThanOnePoint),meanValuesForBinsTET(moreThanOnePoint),stdErrValuesForBinsTET(moreThanOnePoint),...
    '-','Color',MYCOLORS(1,:),'MarkerSize',5,'LineWidth',2);
moreThanOnePoint=countsTemp>1;
errorbar(binCentersTemp(moreThanOnePoint),meanValuesForBinsTemp(moreThanOnePoint),stdErrValuesForBinsTemp(moreThanOnePoint),...
    '-','Color',MYCOLORS(2,:),'MarkerSize',5,'LineWidth',2);
moreThanOnePoint=countsSUL>1;
errorbar(binCentersSUL(moreThanOnePoint),meanValuesForBinsSUL(moreThanOnePoint),stdErrValuesForBinsSUL(moreThanOnePoint),...
    '-','Color',MYCOLORS(3,:),'MarkerSize',5,'LineWidth',2);
if exist('ALSOPLOTDELTAMIN','var')
    moreThanOnePoint=countsDMin>1;
    errorbar(binCentersDMin(moreThanOnePoint),meanValuesForBinsDMin(moreThanOnePoint),stdErrValuesForBinsDMin(moreThanOnePoint),...
        '-','Color',MYCOLORS(4,:),'MarkerSize',5,'LineWidth',2);
end

% Perform & plot tests
% ===

% Perform t-testing
% Note that bins are equal
% Comparing TET with Temp
if ~exist('HIDESQUARES','var')
    for i=find(enoughPointsTET&enoughPointsTemp)
        outcomeNotSameTemp=ttest2(dataPerBinsTET{i},dataPerBinsTemp{i});
        if outcomeNotSameTemp
            plot(binCentersTemp(i),meanValuesForBinsTemp(i),'s','Color',MYCOLORS(2,:),'MarkerSize',15);%,'MarkerFaceColor',MYCOLORS(2,:));
        end
    end
    % Comparing TET with Sul
    for i=find(enoughPointsTET&enoughPointsSUL)
        outcomeNotSameSUL=ttest2(dataPerBinsTET{i},dataPerBinsSUL{i});
        if outcomeNotSameSUL
            plot(binCentersSUL(i),meanValuesForBinsSUL(i),'s','Color',MYCOLORS(3,:),'MarkerSize',15);%,'MarkerFaceColor',MYCOLORS(3,:));
        end
    end
    % Comparing TET with DMin
    if exist('ALSOPLOTDELTAMIN','var')
    for i=find(enoughPointsTET&enoughPointsDMin)
        outcomeNotSameDMin=ttest2(dataPerBinsTET{i},dataPerBinsDMin{i});
        if outcomeNotSameDMin
            plot(binCentersDMin(i),meanValuesForBinsDMin(i),'s','Color',MYCOLORS(4,:),'MarkerSize',15);%,'MarkerFaceColor',MYCOLORS(3,:));
        end
    end
    end
end

% Cosmetics

ylim([0,100]);
xlabel('Birth size (µm)');
ylabel('Interdivision time (min)');
MW_makeplotlookbetter(20);

if ~exist('ALSOPLOTDELTAMIN','var')
    legend({'Tetracycline', 'Temperature', 'SulA'});
else
    legend({'Tetracycline', 'Temperature', 'SulA','Min mutant'});
end
legend boxoff  

hFigSIIC = h;


%% Compare donachie deltamin w. our deltamin

myColors=linspecer(5);

hInterdivMin=figure(6); clf; hold on;
scatter(combinedDynamicsData.('deltaMinTET').selectedXdata,combinedDynamicsData.('deltaMinTET').selectedYdata,4^2,[.7 .7 .7],'filled','MarkerFaceAlpha',.5);
errorbar(combinedDynamicsData.('deltaMinTET').binCentersDynData,combinedDynamicsData.('deltaMinTET').meanValuesForBinsDynData,...
    combinedDynamicsData.('deltaMinTET').stdErrValuesForBins,'o','LineWidth',2,'Color',myColors(4,:),'MarkerSize',7);
plot(donachieInterdivData(:,1),donachieInterdivData(:,2),'o','LineWidth',2','Color',myColors(5,:),'MarkerSize',10);
xlabel('Birth length');
ylabel('Interdivision time');
MW_makeplotlookbetter(20);
title('Min mutant strain');
legend({'This work','This work','Donachie & Begg'});

