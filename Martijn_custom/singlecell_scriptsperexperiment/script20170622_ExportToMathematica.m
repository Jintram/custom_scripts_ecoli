

%%

% Run 
% script20161018_plotAllCrossCorrs_CRP
% first. Specifically the section "Actually plotting the cross-corrs Daan
% style: do it using loop" and the section "".

% Note that the field indices are hard-coded below, so make sure this is
% correct.
%% Convert to output for Laurens' sheet

for groupIdx= 1:numel(GROUPSTOPLOT)
    currentGroupName = GROUPSTOPLOT{groupIdx};
    
    % First cross-correlations
    % ===
    
    % cross-correlation functions for CRP
    CC_conc_growth_CRP    = [gatheredCCs.(currentGroupName).data.Y.datatau{1}; gatheredCCs.(currentGroupName).data.Y.datacorrelation{1}]';        
    CC_prod_growth_CRP    = [gatheredCCs.(currentGroupName).data.Y.datatau{2}; gatheredCCs.(currentGroupName).data.Y.datacorrelation{2}]';
    % cross-correlation functions for constitutive reporter
    CC_conc_growth_consti = [gatheredCCs.(currentGroupName).data.C.datatau{1}; gatheredCCs.(currentGroupName).data.C.datacorrelation{1}]';        
    CC_prod_growth_consti = [gatheredCCs.(currentGroupName).data.C.datatau{2}; gatheredCCs.(currentGroupName).data.C.datacorrelation{2}]';
        
    % cross-correlation function for CRP-const.
    CY_conc_consti_CRP    = [gatheredCYs.(currentGroupName).data.CY.datatau{1}; gatheredCYs.(currentGroupName).data.CY.datacorrelation{1}]';
    CY_prod_consti_CRP    = [gatheredCYs.(currentGroupName).data.CY.datatau{2}; gatheredCYs.(currentGroupName).data.CY.datacorrelation{2}]';
    
    % coss-correlation function for prod-concentration within same fluor
    CC_prodConcConsti = [gatheredPEs.(currentGroupName).data.C.datatau{1}; gatheredPEs.(currentGroupName).data.C.datacorrelation{1}]';
    CC_prodConcCRP    = [gatheredPEs.(currentGroupName).data.Y.datatau{1}; gatheredPEs.(currentGroupName).data.Y.datacorrelation{1}]';
    
    % autocorrelations        
    % ===
    AC_conc_CRP     = [gatheredACs.(currentGroupName).data.Y.datatau{1}; gatheredACs.(currentGroupName).data.Y.datacorrelation{1}]';        
    AC_prod_CRP     = [gatheredACs.(currentGroupName).data.Y.datatau{2}; gatheredACs.(currentGroupName).data.Y.datacorrelation{2}]';
    AC_conc_consti  = [gatheredACs.(currentGroupName).data.C.datatau{1}; gatheredACs.(currentGroupName).data.C.datacorrelation{1}]';        
    AC_prod_consti  = [gatheredACs.(currentGroupName).data.C.datatau{2}; gatheredACs.(currentGroupName).data.C.datacorrelation{2}]';
    AC_growth       = [gatheredACs.(currentGroupName).data.Y.datatau{3}; gatheredACs.(currentGroupName).data.Y.datacorrelation{3}]';
        
    % growth rates
    % ===
    mean_growth_rate = 1/gatheredCYs.(currentGroupName).data.CY.hrsPerDoublingMean;
    
    % save those correlations
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'CC_conc_growth_CRP.mat'],    'CC_conc_growth_CRP');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'CC_prod_growth_CRP.mat'],    'CC_prod_growth_CRP');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'CC_conc_growth_consti.mat'], 'CC_conc_growth_consti');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'CC_prod_growth_consti.mat'], 'CC_prod_growth_consti');
    
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'CC_prodConcConsti.mat'],    'CC_prodConcConsti');    
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'CC_prodConcCRP.mat'],       'CC_prodConcCRP');
    
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'AC_conc_CRP.mat'],     'AC_conc_CRP');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'AC_prod_CRP.mat'],     'AC_prod_CRP');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'AC_conc_consti.mat'],  'AC_conc_consti');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'AC_prod_consti.mat'],  'AC_prod_consti');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'AC_growth.mat'],       'AC_growth');
    
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'CY_conc_consti_CRP.mat'], 'CY_conc_consti_CRP');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'CY_prod_consti_CRP.mat'], 'CY_prod_consti_CRP');
    
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_' 'mean_growth_rate_dbl_per_hr.mat'],  'mean_growth_rate');
        
    % put in all growth rates    
    %{
    paramOfInterestIdx=1; % i.e. growth rate
    nrDatasets=numel(gatheredOutput.(currentGroupName).allData); 
    allGrowthRates=[]; allTimes=[];
    for datasetIdx=1:nrDatasets
        allGrowthRates = [allGrowthRates, ...
            gatheredOutput.(currentGroupName).allData{datasetIdx}{paramOfInterestIdx}'];
        allTimes = [allTimes, ...
            gatheredOutput.(currentGroupName).allData{datasetIdx}{paramOfInterestIdx}'];
    end
    allGrowthRates = allGrowthRates(~isnan(allGrowthRates));
    m3 = allGrowthRates;
    %}
    
end

winopen('U:\Mathematica\CRP_LaurensKrah\data\');
disp('Laurens export mathematica script done..');
