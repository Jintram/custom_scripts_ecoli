
% Make a script that shows the branches and determines the recovery
% points..

%% First manual calculation of experimental switch times

% dataset1 (2013-12-09, pos3)
% G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-09
% first frame   (532): 2013:12:09 17:39:56
% switch frame (1169): 2013:12:10 08:30:54
myTimes(1)=24*60-(17*60+39+56/100) + 08*60+30+54/100
% OR DOES HE MEAN SWITCHTIME = 1169 MINUTES??

% dataset2 (2013-09-24, pos4)
% G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24
% first frame  (728): 2013:09:24 16:48:36
% switch frame (996): 2013:09:24 23:33:11
myTimes(2)=23*60+33+11/100 - (16*60+48+36/100)

% dataset3 (\2013-12-16\pos4crop)
% G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16
% first frame (898): 2013:12:17 15:31
% switch time      : 17-12-2013 15.30 % from exp. description
myTimes(3) = 0

% dataset4 (2013-09-24\pos5crop)
% G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24
% first frame (656): 2013:09:24 14:43:35
% switch frame (996): 2013:09:24 23:33:11
myTimes(4) = 23*60+33+11/100 - (14*60+43+35/100)

% dataset5 (2013-12-16\pos5crop)
% G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16
% first frame 898 (pos5crop): 2013:12:17 15:32:44
% switch time           : 17-12-2013 15.30 % from exp. description
myTimes(5) = 0


%%
OUTPUTFOLDER = 'D:\Local_Data\Dropbox\Dropbox\Filamentation recovery\MW\figures_new\matlab_export\';
SUBFOLDER = 'rutgerdatabranches\';
if ~exist([OUTPUTFOLDER SUBFOLDER],'dir'), mkdir([OUTPUTFOLDER SUBFOLDER]), end

if ~exist('ONESTOANALYZE','var')
    ONESTOANALYZE=[1:11];
end
ONESTOPLOT=[1:5];
datasetsPaths = ...
{ ...
... Note that the parameter ONESTOANALYZE makes a subselection of this data.
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos3_long.mat',...
        ... ^ Comes from ..\2013-12-09\pos3crop\data\
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos4.mat',...
        ... ^ Comes from ..\2013-09-24\pos4crop\data\ (ugly images)
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos4_long.mat',...
        ... ^ Comes from ..\2013-12-16\pos4crop\data\ (ugly images)
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos5.mat',...
        ... ^ Comes from ..\2013-09-24\pos5crop\data\ (very little data)
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos5_long.mat',...
        ... ^ Comes from ..\2013-12-16\pos5crop\data\ (ugly images)
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\2uM_pos2.mat',...
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\2uM_pos4.mat',...
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\2uM_pos6.mat',...
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\10uM_pos1.mat',...
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\10uM_pos3.mat',...
'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\10uM_pos6_long.mat',...        
};
datasetsPaths={datasetsPaths{ONESTOANALYZE}}

datasetNames=...
    {'1uM_pos3_long.mat','1uM_pos4.mat','1uM_pos4_long.mat','1uM_pos5.mat',...
    '1uM_pos5_long.mat',...
    '2uM_pos2.mat','2uM_pos4.mat','2uM_pos6.mat',...
    '10uM_pos1.mat','10uM_pos3.mat','10uM_pos6_long.mat'
    }; 
datasetNames={datasetNames{ONESTOANALYZE}};
%% 

gatheredOutput = {};
for ii = 1:numel(datasetsPaths) %(datasetsPaths)
    
    disp(['plotting schnitz # ' num2str(ii)]);    

    schnitzcells = loadandrename(datasetsPaths{ii});
    
    p.fitTime=[0,2000];
    avgLineColor='r';
    whichParamsToPlot= {{'time', 'muP15_fitNew_all'}};
    fluorColorToPlot=''; theYLabels={'Growth rate [dbl/hr]'};
    
    [h1, gatheredBranchOutput]=MW_plotting_branches_mystyle(p,schnitzcells,avgLineColor,whichParamsToPlot,fluorColorToPlot,theYLabels);
        % outputs figure handle

    title(datasetNames{ii},'Interpreter','None');
        
    saveas(h1, [OUTPUTFOLDER SUBFOLDER 'SVG_' 'branches_dataset_' num2str(ii) '.svg']);
    saveas(h1, [OUTPUTFOLDER SUBFOLDER 'FIG_' 'branches_dataset_' num2str(ii) '.fig']);
    saveas(h1, [OUTPUTFOLDER SUBFOLDER 'PNG_' 'branches_dataset_' num2str(ii) '.png']);
    
    % save for later
    gatheredOutput{ii} = gatheredBranchOutput(1);
    
end


%% determine WT growth rate from 9th dataset
schnitzcells = loadandrename(datasetsPaths{9});

lastSchnitzesLastTimeIndex = find([schnitzcells.time]==max([schnitzcells.time]));
growthData = [schnitzcells.('muP15_fitNew_all')];
lastGrowthRates = growthData(lastSchnitzesLastTimeIndex);
meanGrowthRateAfterRecovery = mean(lastGrowthRates(~isnan(lastGrowthRates)))

% figure(24); xlim([0,2000]);

%% Now find points at which data shows low value last time
%targetValue = meanGrowthRateAfterRecovery/2;
targetValue = meanGrowthRateAfterRecovery/10;

figure(); clf; 
subplot(2,1,1); hold on;

switchPoints={};
for ii=1:numel(datasetsPaths) % ONESTOPLOT
    plot(gatheredOutput{ii}.binCenters,gatheredOutput{ii}.meanValuesForBins,'LineWidth',2)

    %higherValuesIdx     = gatheredOutput{ii}.meanValuesForBins>targetValue;
    lowerValuesIdx      = gatheredOutput{ii}.meanValuesForBins<targetValue;

    %targetIdx = find( (higherValuesIdx(2:end) & lowerValuesIdx(1:end-1)) );
    allTargetIdx = find( lowerValuesIdx );
    targetIdx    = allTargetIdx(end);

    hitTimes  = gatheredOutput{ii}.binCenters(targetIdx);
    hitValues = gatheredOutput{ii}.meanValuesForBins(targetIdx);
    
    plot(hitTimes, hitValues,'o','LineWidth',2);    
    
    switchPoints{ii} = [hitTimes; hitValues];
end

ylim([0,1.5]);

xlabel('time, t [s]');
ylabel('growth rate [dbl/hr]');
MW_makeplotlookbetter(15);


%% Now re-align the graphs on the last times

%figure(); clf; hold on;
subplot(2,1,2); hold on;

shiftTimes=[];
for ii=1:numel(datasetsPaths) %(datasetsPaths)
    
    shiftTimes(ii) = switchPoints{ii}(1,end);    
    plot(gatheredOutput{ii}.binCenters-shiftTimes(ii),gatheredOutput{ii}.meanValuesForBins,'LineWidth',2);
    
end

ylim([0,1.5]);

xlabel('time, t [min]');
ylabel('growth rate [dbl/hr]');
MW_makeplotlookbetter(15);

shiftTimes

%%
hTraces= figure; clf; hold on;

shiftTimes=[]; % could use "myTimes" instead to get actual switch times
for ii=1:numel(datasetsPaths) %(datasetsPaths)
    
    shiftTimes(ii) = switchPoints{ii}(1,end);    
    plot(gatheredOutput{ii}.binCenters-myTimes(ii),gatheredOutput{ii}.meanValuesForBins,'LineWidth',2);
    
end

ylim([0,1.1]);
xlim([-1250,500]);

xlabel('Time (min)');
ylabel('Growth rate (dbl/hr)');
MW_makeplotlookbetter(20);