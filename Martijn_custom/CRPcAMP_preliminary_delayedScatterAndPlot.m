
%% Description

% MW 2015/04
%
% Script to quickly generate "delayed scatter" plots, for now based on
% NW's data handling; 
% MW TODO: Check PN data.
%
% About input params:
% - associatedFieldNames{1} MUST BE timefield field name
% - associatedFieldNames{2} X, usually fluor concentration or rate field name
% - associatedFieldNames{3} Y, usually growth rate field name
%
% - makeDtime   if set to 0 nothing happens, if set to 1, extra 
%               MWDJK_dX_time filed is recalculated.
%
% Example of how to call script:
%{
myID = ' icd-lac';
p.movieName = 'pos6crop';
p.movieDate = '2013-05-16';

myOutputFolder = ['F:\X_Other_datasets\CRPcAMP\NWdata-pfkA\' myID '_' p. movieDate  '_' p.movieName '\'];

myDataFile = 'F:\X_Other_datasets\CRPcAMP\NWdata-pfkA\icd-lac-2013-05-16-pos6crop-Schnitz.mat';
associatedFieldNames =  {'G_time','G6_mean', 'muP15_fitNew'};
myTitle = ['NW' myID '_' p. movieDate  '_' p.movieName];

myFitTime = [50   1000];

p.NW_saveDir = [myOutputFolder 'misc\'];
p.DJK_saveDir = [myOutputFolder 'misc\'];

p.fluor1='r';
p.fluor2='g';

% info required to make branches
badSchnitzes = [428]; % from excel top field “schnitzes to be removed from analysis“

CRPcAMP_preliminary_delayedScatterAndPlot

   
associatedFieldNames =  {'dG5_time','dG5_cycCor', 'muP15_fitNew'};
CRPcAMP_preliminary_delayedScatterAndPlot
%}


if ~exist('myOutputFolder')
    myOutputFolder = 'C:\Users\wehrens\Desktop\testdelayedscatter\output\';
end

if ~exist('associatedFieldNames') | ~exist('myTitle') | ~exist('p') | ~exist('badSchnitzes')
    error('input not supplied.')
end
    

% Loading
%load(myDataFile); % Can also be done by user


%%

name_rm = 'rm'; name_all = 'all';
fitTime = myFitTime;

if ~alreadyRemovedInMatFile
    % Find Schnitzes with slow/negative growth rate -> rm them?!
    slowschnitzes=NW_detectSlowSchnitzes(p,schnitzcells,associatedFieldNames{3},'muThreshold',0.05);

    % Preparation to load
    % Adapted from NW excel sheet    
    s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); 
    s_all_fitTime = DJK_selSchitzesToPlot(s_all, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); name_all_fitTime = ['all_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
    s_all_fitTime_cycle = DJK_selSchitzesToPlot(s_all_fitTime, 'completeCycle', @(x) x ~= 0); name_all_fitTime_cycle = [name_all_fitTime '_cycle'];

    s_rm = DJK_selSchitzesToPlot(s_all, 'P', @(x) 1); 
    for i=badSchnitzes, s_rm(i).useForPlot=0; end;
    s_rm_fitTime = DJK_selSchitzesToPlot(s_rm, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); name_rm_fitTime = ['rm_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
    s_rm_fitTime_cycle = DJK_selSchitzesToPlot(s_rm_fitTime, 'completeCycle', @(x) x ~= 0); name_rm_fitTime_cycle = [name_rm_fitTime '_cycle'];
end

%% 
% Collect indices of the GFP measurement per schnitz.
% This is not necessary because Noreen's data already contains the right
% fields -and the fluor code already does this 
% (DJK_addToSchnitzes_fluor_anycolor) - but might be convenient for future 
% generation of "fieldX_at_G".
indicesForFluor = {};
FluorFieldAllName = [upper(p.fluor1) '_mean_all'];
FluorRateBaseFieldAllName = [upper(p.fluor1) '_time'];
numelS_rm = numel(s_rm);
for i = 1:numelS_rm
    % current schnitz
    s = s_rm(i);
    
    % Re-calculate fluor times
    indicesForFluor{i} = find(~isnan(s_rm(i).(FluorFieldAllName))); % could access field name as string
    % This is only valid for concentration values:
    theTimes = s_rm(i).('time');
    s_rm(i).MWTimeField = theTimes(indicesForFluor{i});
    
    % Calculate dX_time field if it doesn't exist, based on the assumption
    % it was made by DJK_addTocSchnitzes_fluorRate_phase, which takes into
    % acount the points t+1 and t.
    if makeDtime        
            
            % include the daughter field that was also used for fluor calc
            if s.D ~= 0
                % take only one daughter for time reference
                sD = s_rm(s.D);
                % extend
                extendedTime = [s.(FluorRateBaseFieldAllName), sD.(FluorRateBaseFieldAllName)(1)];
            else
                % end of lineage, no extension
                extendedTime = [s.(FluorRateBaseFieldAllName)]; 
            end
            
            % calculate times
            % sum n and n+1, average
            DTimes = [extendedTime(1:end-1)+extendedTime(2:end)]./2;
                    
            % add to schnitz
            fieldName = ['MWDJK_' 'd' FluorRateBaseFieldAllName];
            s_rm(i).(fieldName) = DTimes;

    end
end

if makeDtime
    warning('Determination of dX_time values only valid for DJK_addToSchnitzes_fluor_anycolor rates.');
    pause(1);
end

%% Calculating branches
% ===
s_rm = MW_calculateframe_nrs(s_rm); % backwards compatibility fix 

fitTime = fitTime + [2 -2];


branchData = DJK_getBranches(p,s_rm,'dataFields',{associatedFieldNames{1}, associatedFieldNames{2}, associatedFieldNames{3} }, 'fitTime', fitTime); 
 name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2)) '_Conc_oldRates'];


%%

% Just some plot colors
distinguishableColors = distinguishable_colors(numel(branchData)+1,[1 1 1]); 

% Plot all branches
figure(1); clf; hold on;
numelBranches = numel(branchData);
for branch_nr = 1:numelBranches
    l = plot(branchData(branch_nr).(associatedFieldNames{1}), branchData(branch_nr).(associatedFieldNames{3}),'-o','Color',distinguishableColors(branch_nr,:))
    set(l, 'LineWidth', (numelBranches-branch_nr+1)/numelBranches*10);
end

xlabel(associatedFieldNames{1},'Interpreter', 'None'), ylabel(associatedFieldNames{3},'Interpreter', 'None')

%%

REDUNDANCYALLOWED = 2^2;

% Some additional editing of the branches:
branchData = DJK_addToBranches_noise(p, branchData,'dataFields',{associatedFieldNames{1},associatedFieldNames{2},associatedFieldNames{3}});
%trimmed_branches = DJK_trim_branch_data(branches);
branch_groups = DJK_divide_branch_data(branchData);

% Colony average mean has already been substracted so theoretically extra
% normalization shouldn't have an effect.
p.extraNorm=0;

% THIS MIGHT FAIL BECAUSE TIME IS NOT SET CORRECTLY (SHOULD BE time_at_..)
% To calculate cross correlations, additional normalization is usually
% performed, namely to filter out colony average behavior.
[CorrData,composite_corr] = DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, ['noise_' associatedFieldNames{1,2}],['noise_' associatedFieldNames{1,3}] ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',1); 

% Do we want to filter out colony average behavior for the "delayed
% scatter" plots also? Maybe do this with noise fields?
% But let's try with "raw" data first..
p.timeField = associatedFieldNames{1,1};
% p.tauIndices = [-7:7]; %p.tauIndices = [-29:4:-1,0,1:4:30];
if isfield(p,'tauIndices'), p=rmfield(p,'tauIndices'); end
[dataPairsPerTau, iTausCalculated, originColorPerTau, correlationsPerTau] = ...
    MW_getdelayedscatter(p, branchData, ['noise_' associatedFieldNames{1,2}], ['noise_' associatedFieldNames{1,3}], REDUNDANCYALLOWED)

%% Compare two cross-corrs (DJK & MW)

figure(1),clf,hold on
errorbar(CorrData(:,1),CorrData(:,2),CorrData(:,3),'x-','Color', [.5,.5,.5], 'LineWidth',2)
l1=plot(CorrData(:,1),CorrData(:,2),'x-k','LineWidth',2)

l2=plot(CorrData(:,1),correlationsPerTau,'o-r','LineWidth',2)
myxlimvalues=[min(CorrData(:,1)), max(CorrData(:,1))];
xlim(myxlimvalues);
ylim([-0.4,0.4]);
plot(myxlimvalues,[0,0],'k-');

legend([l1,l2],{'DJK','MW'})

title(['DJK vs. MW R -- ' myID '_' p. movieDate  '_' p.movieName], 'Interpreter', 'none');
xlabel('\tau (hrs)');
ylabel(['R(' associatedFieldNames{1,2} ', growth) (normalized)'], 'Interpreter', 'none');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

saveas(1,[myOutputFolder 'crosscorrs_' associatedFieldNames{1,2} '.tiff']);



%% Plot code from CRPcAMP..overview..general
% ==========
NRCONTOURLINES = 5;
SHOWPLUSMINFROMZERO = 15;
PLOT3DSCATTER = 0;

middlePosition = [ceil(numel(iTausCalculated)/2)-SHOWPLUSMINFROMZERO:ceil(numel(iTausCalculated)/2)+SHOWPLUSMINFROMZERO];
rangeiTausCalculated = middlePosition;
% delayIdx = 11; % 39 is middle

myColorMap = colormap(winter(numel(iTausCalculated)));

if PLOT3DSCATTER
    hFig = figure(2); clf; hold on;
    offset=100; width1=500; height1=500;
    set(hFig, 'Position', [offset offset width1 height1]);
end

hFig = figure(3); clf; hold on;

densities=[]; Xs=[]; Ys=[]; Zs=[];
for delayIdx = rangeiTausCalculated
    
    % Rename data more convenient
    data = [dataPairsPerTau{delayIdx}(:,1), dataPairsPerTau{delayIdx}(:,2)];    
    
    % plot scatter
    if PLOT3DSCATTER
        hFig = figure(2), hold on;
        scatter3(data(:,1),data(:,2),ones(1,numel(data(:,2)))*iTausCalculated(delayIdx),3,myColorMap(delayIdx,:),'.');%,'Color',distinguishableColors(myGrouping(i)+1,:),'MarkerSize',3);
    end

    % plot
    set(0,'CurrentFigure',3); clf, hold on;
    offset=100; width1=500; height1=500;
    set(hFig, 'Position', [(offset+width1) offset width1 height1]);     
     
    % scatter
    for pointIdx = 1:numel(data(:,1))
        plot(data(pointIdx,1),data(pointIdx,2),'.','Color',originColorPerTau{delayIdx}(pointIdx,:));%,'Color',distinguishableColors(myGrouping(i)+1,:),'MarkerSize',3);
    end
    
    % contour (from kde)
    [bandwidth,density,X,Y] = kde2d(data);    
    [C, l1] = contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);
          
    % mean
    lineH = plot(mean(data(:,1)),mean(data(:,2)),'o','MarkerFaceColor','k','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);
    
    title(['D# = ' num2str(iTausCalculated(delayIdx)) ', R = ' sprintf('%0.3f',correlationsPerTau(delayIdx)), ' -- ' myID '_' p. movieDate  '_' p.movieName], 'Interpreter', 'none')    
    
    xlabel(['Delta ' associatedFieldNames{1,2}] , 'Interpreter', 'none');
    ylabel(['Delta ' associatedFieldNames{1,3}] , 'Interpreter', 'none');
    %Set all fontsizes
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
    set(gca,'FontSize',15);
    
    midIdx = ceil(numel(iTausCalculated)/2);
    xlim([  min(dataPairsPerTau{midIdx}(:,1)), max(dataPairsPerTau{midIdx}(:,1))  ])
    ylim([  min(dataPairsPerTau{midIdx}(:,2)), max(dataPairsPerTau{midIdx}(:,2))  ])
       
    saveas(3,[myOutputFolder 'graphTauIdx_' associatedFieldNames{1,2} '_' sprintf('%05d',delayIdx) '.tiff']);
   
end

figure(3);

% average point (used for legend too)
%{
legendLines = []; previous = 0;
for i = 1:numberOfDataFiles        
    lineH = plot(mean(data(:)),mean(data(2,:)),'o','MarkerFaceColor',distinguishableColors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);
if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','northeast');
%}

if PLOT3DSCATTER
    figure(2);
    xlabel('Growth rate (dbl/hr)');
    ylabel('Concentration (a.u.)');
    %Set all fontsizes
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
    set(gca,'FontSize',15);
    %ylim([-750, 2000])
    %xlim([0, max([myData(:).selected_growth_rates])])
end
    

%% Some random checks

% distribution of times 
% equivalent of 'hist' calculated manually because times have interval,
% which complicates use of hist function.
fieldcount=0; R = struct;
% loop over timefields of interest
for myfield = {'Y_time','MWDJK_dY_time', 'time'}

    % administration
    myfield = char(myfield);
    fieldcount = fieldcount+1;
    
    % get all possible values 
    R(fieldcount).YValues = unique([s_rm.(myfield)]);
    % count how many of each values are in the schnitz struct
    R(fieldcount).countYValues = [];
    for value=R(fieldcount).YValues 
        R(fieldcount).countYValues(end+1) = sum([s_rm.(myfield)]==value);
    end
    
end
% plot first two fields of interest
figure, clf; hold on;
if numel(R)>2, plot(R(3).YValues, R(3).countYValues, 'xk'); end
plot(R(1).YValues, R(1).countYValues, 'or','LineWidth',2);
plot(R(2).YValues, R(2).countYValues, 'ob','LineWidth',2);








