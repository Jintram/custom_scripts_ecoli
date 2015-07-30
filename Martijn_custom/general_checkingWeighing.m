% This script can be run in combination with CRPcAMP_preliminary_delayedScatterAndPlot

% Testing whether the weighing procedure can introduce an artificial correlation
somerandomXdata = rand(1000,1);
somerandomYdata = rand(1000,1);

Rnoise = corr(somerandomXdata,somerandomYdata)

% First make our schnitzes have random values for the growth speeds
uncorrelated_schnitz = s_rm;

for i = 1:numel(uncorrelated_schnitz)   
    lengthFields=numel(uncorrelated_schnitz(i).G_time);
    uncorrelated_schnitz(i).whitenoiseX = rand(1,lengthFields); 
    uncorrelated_schnitz(i).whitenoiseY = rand(1,lengthFields);
end

% Generate branches on oncorrelated schnitz data
branchData = DJK_getBranches(p,uncorrelated_schnitz,'dataFields',{'G_time', 'whitenoiseX', 'whitenoiseY' }, 'fitTime', fitTime); 

% Substract colony average
branchData = DJK_addToBranches_noise(p, branchData,'dataFields',{associatedFieldNames{1},'whitenoiseX','whitenoiseY'});

% create branch_groups
branch_groups = DJK_divide_branch_data(branchData);

% THIS MIGHT FAIL BECAUSE TIME IS NOT SET CORRECTLY (SHOULD BE time_at_..)
% To calculate cross correlations, additional normalization is usually
% performed, namely to filter out colony average behavior.
[CorrData,composite_corr] = DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_whitenoiseX','noise_whitenoiseY' ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',1);



