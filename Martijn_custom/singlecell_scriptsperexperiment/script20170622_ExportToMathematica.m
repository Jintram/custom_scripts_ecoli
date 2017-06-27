

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

    % convert to matrices as desired by Laurens
    % first data is R_concentration-mu, hence {1}
    m1 = [gatheredCCs.(currentGroupName).data.datatau{1}; gatheredCCs.(currentGroupName).data.datacorrelation{1}]';    
    % second data is R_production-mu, hence {2}
    m2 = [gatheredCCs.(currentGroupName).data.datatau{2}; gatheredCCs.(currentGroupName).data.datacorrelation{2}]';
    % put in all growth rates    
    paramOfInterestIdx=1; % i.e. growth rate
    nrDatasets=numel(gatheredOutput.(currentGroupName).allData); 
    allGrowthRates=[];
    for datasetIdx=1:nrDatasets
        allGrowthRates = [allGrowthRates, ...
            gatheredOutput.(currentGroupName).allData{datasetIdx}{paramOfInterestIdx}'];
    end
    allGrowthRates = allGrowthRates(~isnan(allGrowthRates));
    m3 = allGrowthRates;
    
    % save those matrices
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_R_concentration_lambda.mat'],    'm1');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_R_production_lambda.mat'],       'm2');
    save(['U:\Mathematica\CRP_LaurensKrah\data\' currentGroupName '_Growth_rates_pile.mat'],         'm3');
    
end
