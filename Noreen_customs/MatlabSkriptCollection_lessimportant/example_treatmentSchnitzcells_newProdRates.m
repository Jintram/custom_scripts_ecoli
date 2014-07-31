close all
clear all
p = DJK_initschnitz('pos1crop','2012-05-18','e.coli.AMOLF','rootDir','D:\colonies\', 'cropLeftTop', [181,105], 'cropRightBottom', [1090,1020],'fluor1','y','fluor2','r','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);
munb = '17';

schnitzcells = PN_fluorRate(schnitzcells);
schnitzcells = PN_fluorRate_R(schnitzcells);

% CYCLE CORRECT FLUOR RATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
schnitzcells = DJK_addToSchnitzes_predictedValues(schnitzcells, 'phase', 'length_fitNew', 'phase2', [0 1]);
schnitzcells = DJK_addToSchnitzes_atYandDY(schnitzcells, 'phase2');
schnitzcells = PN_addToSchnitzes_Phase_at_TimeField(schnitzcells,'phase2','dY5_time');
schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells, 'dY5', 'phase2_at_dY5_time');
schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells, 'dR5', 'phase2_at_dY5_time');
schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells, 'R6_mean', 'phase2_atY');
schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells, 'Y6_mean', 'phase2_atY');
schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells, ['muP' munb '_fitNew'], 'phase2_atY');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fitTime = [400 1300]; 
maxlim = 1300;
autolim = [0 maxlim];
crosslim = [-maxlim maxlim];
rmschnitz = [390 454 456];
fitTime = fitTime + [2 -2];

s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';
s_all_fitTime = DJK_selSchitzesToPlot(s_all, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); name_all_fitTime = ['all_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
s_all_fitTime_cycle = DJK_selSchitzesToPlot(s_all_fitTime, 'completeCycle', @(x) x ~= 0); name_all_fitTime_cycle = [name_all_fitTime '_cycle'];
s_rm = DJK_selSchitzesToPlot(s_all, 'P', @(x) 1); name_rm = 'rm';
for i=[rmschnitz], s_rm(i).useForPlot=0; end;
s_rm_fitTime = DJK_selSchitzesToPlot(s_rm, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); name_rm_fitTime = ['rm_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
s_rm_fitTime_cycle = DJK_selSchitzesToPlot(s_rm_fitTime, 'completeCycle', @(x) x ~= 0); name_rm_fitTime_cycle = [name_rm_fitTime '_cycle'];
% 

datafields = {'dY5_time' 'dR5_cycCor' 'dY5_cycCor' ['muP' munb '_fitNew_cycCor']};
branchData = DJK_get_branches(p,s_rm,'dataFields',datafields, 'fitTime', fitTime); name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
DJK_xcorr_branches(p, branchData, 'noise_dY5_cycCor', 'noise_dY5_cycCor', 'correct', 1, 'xlim', autolim, 'selectionName', name_rm_branch, 'onScreen', 0,'timeField','dY5_time');
DJK_xcorr_branches(p, branchData, 'noise_dY5_cycCor', ['noise_muP' munb '_fitNew_cycCor'], 'correct', 1, 'xlim',crosslim, 'selectionName', name_rm_branch, 'onScreen', 0,'timeField','dY5_time');
DJK_xcorr_branches(p, branchData, 'noise_dR5_cycCor', 'noise_dR5_cycCor', 'correct', 1, 'xlim', autolim, 'selectionName', name_rm_branch, 'onScreen', 0,'timeField','dY5_time');
DJK_xcorr_branches(p, branchData, 'noise_dR5_cycCor', ['noise_muP' munb '_fitNew_cycCor'], 'correct', 1, 'xlim',crosslim, 'selectionName', name_rm_branch, 'onScreen', 0,'timeField','dY5_time');
DJK_xcorr_branches(p, branchData, 'noise_dY5_cycCor','noise_dR5_cycCor', 'correct', 1, 'xlim',crosslim, 'selectionName', name_rm_branch, 'onScreen', 0,'timeField','dY5_time');
branchData = DJK_getBranches(p,s_rm,'dataFields',datafields, 'fitTime', fitTime); name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
branches = DJK_addToBranches_noise(p, branchData,'dataFields',datafields);
trimmed_branches = DJK_trim_branch_data(branches);
branch_groups = DJK_divide_branch_data(trimmed_branches);
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dY5_cycCor', 'noise_dY5_cycCor','selectionName',name_rm_branch,'timeField','dY5_time');
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dY5_cycCor', ['noise_muP' munb '_fitNew_cycCor'],'selectionName',name_rm_branch,'timeField','dY5_time');
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dR5_cycCor', 'noise_dR5_cycCor','selectionName',name_rm_branch,'timeField','dY5_time');
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dR5_cycCor', ['noise_muP' munb '_fitNew_cycCor'],'selectionName',name_rm_branch,'timeField','dY5_time');
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dY5_cycCor', 'noise_dR5_cycCor','selectionName',name_rm_branch,'timeField','dY5_time');


datafields = {'Y_time' 'Y6_mean_cycCor' 'R6_mean_cycCor' ['muP' munb '_fitNew_cycCor']};
branchData = DJK_get_branches(p,s_rm,'dataFields',datafields, 'fitTime', fitTime); name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
DJK_xcorr_branches(p, branchData, 'noise_Y6_mean_cycCor', 'noise_Y6_mean_cycCor', 'correct', 1, 'xlim', autolim, 'selectionName', name_rm_branch, 'onScreen', 0);
DJK_xcorr_branches(p, branchData, 'noise_Y6_mean_cycCor', ['noise_muP' munb '_fitNew_cycCor'], 'correct', 1, 'xlim',crosslim, 'selectionName', name_rm_branch, 'onScreen', 0);
DJK_xcorr_branches(p, branchData, 'noise_R6_mean_cycCor', 'noise_R6_mean_cycCor', 'correct', 1, 'xlim', autolim, 'selectionName', name_rm_branch, 'onScreen', 0);
DJK_xcorr_branches(p, branchData, 'noise_R6_mean_cycCor', ['noise_muP' munb '_fitNew_cycCor'], 'correct', 1, 'xlim',crosslim, 'selectionName', name_rm_branch, 'onScreen', 0);
DJK_xcorr_branches(p, branchData, 'noise_Y6_mean_cycCor','noise_R6_mean_cycCor', 'correct', 1, 'xlim',crosslim, 'selectionName', name_rm_branch, 'onScreen', 0);
branchData = DJK_getBranches(p,s_rm,'dataFields',datafields, 'fitTime', fitTime); name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
branches = DJK_addToBranches_noise(p, branchData,'dataFields',datafields);
trimmed_branches = DJK_trim_branch_data(branches);
branch_groups = DJK_divide_branch_data(trimmed_branches);
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_Y6_mean_cycCor', 'noise_Y6_mean_cycCor','selectionName',name_rm_branch);
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_Y6_mean_cycCor',  ['noise_muP' munb '_fitNew_cycCor'],'selectionName',name_rm_branch);
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_R6_mean_cycCor', 'noise_R6_mean_cycCor','selectionName',name_rm_branch);
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_R6_mean_cycCor',['noise_muP' munb '_fitNew_cycCor'],'selectionName',name_rm_branch);
PN_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_Y6_mean_cycCor', 'noise_R6_mean_cycCor','selectionName',name_rm_branch);


% 
% DJK_plot_scatterColor(p, s_all, 'av_mu_fitNew', 'av_time', 'gen', 'ylim', [0 1.0], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);
% DJK_plot_scatterColor(p, s_all, 'av_Y6_mean', 'av_time', 'gen', 'ylim', [0 500], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);
% DJK_plot_scatterColor(p, s_all, 'av_R6_mean', 'av_time', 'gen', 'ylim', [0 300], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);
% DJK_plot_scatterColor(p, s_all, 'av_mu_fitNew', 'av_Y6_mean', 'av_time', 'xlim', [0 500], 'ylim', [0 1.0], 'selectionName', name_all, 'plotRegression', 1, 'onScreen', 0);
% DJK_plot_scatterColor(p, s_all, 'av_mu_fitNew', 'av_R6_mean', 'av_time', 'xlim', [0 300], 'ylim', [0 1.0], 'selectionName', name_all, 'plotRegression', 1, 'onScreen', 0);
% 
% DJK_plot_scatterColor(p, s_all_fitTime_cycle, 'av_mu_fitNew', 'av_Y6_mean', 'av_time', 'xlim', [0 500], 'ylim', [0 1.0], 'selectionName', name_all_fitTime_cycle, 'plotRegression', 1, 'onScreen', 0);
% DJK_plot_scatterColor(p, s_rm_fitTime_cycle, 'av_mu_fitNew', 'av_Y6_mean', 'av_time', 'xlim', [0 500], 'ylim', [0 1.0], 'selectionName', name_rm_fitTime_cycle, 'plotRegression', 1, 'onScreen', 0);
% DJK_plot_scatterColor(p, s_all_fitTime_cycle, 'av_mu_fitNew', 'av_R6_mean', 'av_time', 'xlim', [0 300], 'ylim', [0 1.0], 'selectionName', name_all_fitTime_cycle, 'plotRegression', 1, 'onScreen', 0);
% DJK_plot_scatterColor(p, s_rm_fitTime_cycle, 'av_mu_fitNew', 'av_R6_mean', 'av_time', 'xlim', [0 300], 'ylim', [0 1.0], 'selectionName', name_rm_fitTime_cycle, 'plotRegression', 1, 'onScreen', 0);
% 
% % Time point data
% schnitzData = DJK_get_schnitzData(p, s_rm_fitTime,'dY5_time', 'dataFields',datafields, 'fitTime', fitTime);
% 
% DJK_plot_scatterColor(p, schnitzData,['muP' munb '_fitNew_cycCor'], 'dY5_cycCor', 'dY5_time', 'xlim', [0 300], 'ylim', [0 1.0], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% DJK_plot_scatterColor(p, schnitzData,['muP' munb '_fitNew_cycCor'], 'noise_dY5_cycCor', 'dY5_time', 'xlim', [-50 50], 'ylim', [-0.7 0.7], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% DJK_plot_scatterColor(p, schnitzData,['muP' munb '_fitNew_cycCor'], 'dR5_cycCor', 'dY5_time', 'xlim', [0 300], 'ylim', [0 1.0], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% DJK_plot_scatterColor(p, schnitzData,['muP' munb '_fitNew_cycCor'], 'noise_dR5_cycCor', 'dY5_time', 'xlim', [-50 50], 'ylim', [-0.7 0.7], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% DJK_plot_scatterColor(p, schnitzData, 'dY5_cycCor', 'dR5_cycCor', 'dY5_time', 'xlim', [0 500], 'ylim', [0 500], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% DJK_plot_scatterColor(p, schnitzData, 'noise_dY5_cycCor', 'noise_dR5_cycCor', 'dY5_time', 'xlim', [-200 200], 'ylim', [-500 500], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% 
% DJK_plot_time_hist(p, schnitzData, ['muP' munb '_fitNew_cycCor'], 0, 'binCenters', [0:0.05:1.0], 'selectionName', name_rm_fitTime, 'onScreen', 0,'timeField','dY5_time');
% DJK_plot_time_hist(p, schnitzData, 'dY5_cycCor', 0.94, 'binCenters', [0:25:2000], 'selectionName', name_rm_fitTime, 'onScreen', 0,'timeField','dY5_time');
% DJK_plot_time_hist(p, schnitzData, 'dR5_cycCor', 0.94, 'binCenters', [0:25:2000], 'selectionName', name_rm_fitTime, 'onScreen', 0,'timeField','dY5_time');
% 
% 
% %DJK_plot_scatterColor(p, schnitzData, ['muP' munb '_fitNew_cycCor'], 'Y6_mean_cycCor', 'Y_time', 'xlim', [0 5000], 'ylim', [0 1.0], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% %DJK_plot_scatterColor(p, schnitzData, ['muP' munb '_fitNew_cycCor'], 'noise_Y6_mean_cycCor', 'Y_time', 'xlim', [-50 50], 'ylim', [-0.7 0.7], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% %DJK_plot_time_hist(p, schnitzData, 'Y6_mean_cycCor', 0.94, 'binCenters', [0:5:500], 'selectionName', name_rm_fitTime, 'onScreen', 0);
% %DJK_plot_scatterColor(p, schnitzData,['muP' munb '_fitNew_cycCor'], 'R6_mean_cycCor', 'Y_time', 'xlim', [0 300], 'ylim', [0 1.0], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% %DJK_plot_scatterColor(p, schnitzData,['muP' munb '_fitNew_cycCor'], 'noise_R6_mean_cycCor', 'Y_time', 'xlim', [-50 50], 'ylim', [-0.7 0.7], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% %DJK_plot_time_hist(p, schnitzData, 'R6_mean_cycCor',0,'binCenters', [0:5:300], 'selectionName', name_rm_fitTime, 'onScreen', 0,'timeField','Y_time');
% %DJK_plot_scatterColor(p, schnitzData, 'Y6_mean_cycCor', 'R6_mean_cycCor', 'Y_time', 'xlim', [0 500], 'ylim', [0 500], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
% %DJK_plot_scatterColor(p, schnitzData, 'noise_Y6_mean_cycCor', 'noise_R6_mean_cycCor', 'Y_time', 'xlim', [-200 200], 'ylim', [-500 500], 'selectionName', name_rm_fitTime, 'plotRegression', 1, 'onScreen', 0);
