
%% Description

% MW 2015/04

% Quick script to otbain branchdata/crosscorrs for checking purposes.
%
% (Not a key script.)


%%
% data conversion/loading/renaming
% Section obtained from ecxel file
% "pfkA-lac-schnitzcell_fullAnalysis_2013-04-16pos3.xlsx"
p = DJK_initschnitz('pos3crop','2013-04-16','e.coli.AMOLF','rootDir','\\BIOFYSICASRV\Users2\Walker\Experiments\GrowthNoise\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','r','fluor2','g','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

%slowschnitzes=NW_detectSlowSchnitzes(p,schnitzcells,'muP15_fitNew','muThreshold',0.05);

fitTime = [100 630];
s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';
s_all_fitTime = DJK_selSchitzesToPlot(s_all, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); name_all_fitTime = ['all_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
s_all_fitTime_cycle = DJK_selSchitzesToPlot(s_all_fitTime, 'completeCycle', @(x) x ~= 0); name_all_fitTime_cycle = [name_all_fitTime '_cycle'];


s_rm = DJK_selSchitzesToPlot(s_all, 'P', @(x) 1); name_rm = 'rm';
for i=[ 461 397 485 392 785 792 736 483 567 437 802 618 834 869 862 375 881 895 482 541 546 550 638 683 719], s_rm(i).useForPlot=0; end;
s_rm_fitTime = DJK_selSchitzesToPlot(s_rm, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); name_rm_fitTime = ['rm_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
s_rm_fitTime_cycle = DJK_selSchitzesToPlot(s_rm_fitTime, 'completeCycle', @(x) x ~= 0); name_rm_fitTime_cycle = [name_rm_fitTime '_cycle'];




%% Calculating branches
% ===
s_rm = MW_calculateframe_nrs(s_rm); % backwards compatibility fix 

fitTime = [100 630]; fitTime = fitTime + [2 -2];

% from noreen's excel file
branchData1 = DJK_getBranches(p,s_rm,'dataFields',{'R_time' 'R6_mean' 'G_time','G6_mean' 'muP15_fitNew' }, 'fitTime', fitTime); 
name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2)) '_Conc_oldRates'];
% edited to obtain more fields:
branchData2 = DJK_getBranches(p,s_rm,'dataFields',{'R_time' 'R6_mean' 'G_time','G6_mean' 'muP15_fitNew', associatedFieldNames{1,1},  associatedFieldNames{1,2},  associatedFieldNames{1,3} , 'frame_nrs', 'time'}, 'fitTime', fitTime); 
name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2)) '_Conc_oldRates'];

% Collect indices of the GFP measurement per schnitz.
% This is not necessary because Noreen's data already contains the right
% fields -and the fluor code already does this 
% (DJK_addToSchnitzes_fluor_anycolor) - but might be convenient for future 
% generation of "fieldX_at_G".
indicesForG = {};
for i = 1:numel(s_rm)    
    indicesForG{i} = find(~isnan(s_rm(i).G_mean_all));
end
branchData3 = DJK_getBranches(p,s_rm,'dataFields',{'G_time','G6_mean' 'muP15_fitNew', associatedFieldNames{1,1}, associatedFieldNames{1,2}, associatedFieldNames{1,3}}, 'fitTime', fitTime); 

%%

% Just some plot colors
distinguishableColors = distinguishable_colors(numel(branchData)+1,[1 1 1]); 

% Plot all branches
figure(1); clf; hold on;
numelBranches = numel(branchData);
for branch_nr = 1:numelBranches
    l = plot(branchData(branch_nr).time, branchData(branch_nr).(associatedFieldNames{1}),'-o','Color',distinguishableColors(branch_nr,:))
    set(l, 'LineWidth', (numelBranches-branch_nr+1)/numelBranches*10);
end


%%

% Some additional editing of the branches:
branches = DJK_addToBranches_noise(p, branchData,'dataFields',{'R_time' 'R6_mean' 'G_time' 'G6_mean' 'muP15_fitNew' });
%trimmed_branches = DJK_trim_branch_data(branches);
branch_groups = DJK_divide_branch_data(branches);

% Means are not substracted yet, but the function can be told to do so with
% the following argument:
p.extraNorm=1;

% THIS MIGHT FAIL BECAUSE TIME IS NOT SET CORRECTLY (SHOULD BE time_at_..)
% To calculate cross correlations, additional normalization is usually
% performed, namely to filter out colony average behavior.
DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, associatedFieldNames{1,3},associatedFieldNames{1,1} ,'selectionName',name_rm_branch,'timeField','time','onScreen',1); 

% Do we want to filter out colony average behavior for the "delayed
% scatter" plots also? Maybe do this with noise fields?
% But let's try with "raw" data first..
p.timeField = 'G_time';
[dataPairsPerTau, iTausCalculated] = MW_getdelayedscatter(p, branchData3, associatedFieldNames{1,1}, associatedFieldNames{1,3}, 3)

% Plot data
figure, plot(dataPairsPerTau{39}(:,1)',dataPairsPerTau{39}(:,2)','x')
figure, plot(dataPairsPerTau{30}(:,1)',dataPairsPerTau{30}(:,2)','x')

%% Plot 3d scatter



%% Plot code from CRPcAMP..overview..general
% ==========
NRCONTOURLINES = 2;
delayIdx = 11; % 39 is middle

hFig = figure(2); clf; hold on;
offset=100; width1=500; height1=500;
set(hFig, 'Position', [offset offset width1 height1]);

% Growth rate data
% ===
data = [dataPairsPerTau{delayIdx}(:,1), dataPairsPerTau{delayIdx}(:,2)];
[bandwidth,density,X,Y] = kde2d(data);    
[C, l1] = contour3(X,Y,density,NRCONTOURLINES,'-'); 
set(l1, 'LineWidth', 2,'Color', 'b');
plot(data(:,1),data(:,2),'.');%,'Color',distinguishableColors(myGrouping(i)+1,:),'MarkerSize',3);


% average point (used for legend too)
%{
legendLines = []; previous = 0;
for i = 1:numberOfDataFiles        
    lineH = plot(mean(data(:)),mean(data(2,:)),'o','MarkerFaceColor',distinguishableColors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);
if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','northeast');
%}

xlabel('Growth rate (dbl/hr)');
ylabel('Production rate (a.u./min)');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);
%ylim([-750, 2000])
xlim([0, max([myData(:).selected_growth_rates])])


