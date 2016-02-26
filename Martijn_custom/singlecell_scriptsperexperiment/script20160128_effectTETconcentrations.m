%% Perform analysis

SCRIPTNAME = 'script20160128_effectTETconcentrations';

%{

clear databaseValuesSignalNoise databaseValuesMeanSignalNoise databaseValuesStdSignalNoise databaseValuesNames 

% Data is already saved in .mat files, speficially 
% F:\A_Tans1_step1_incoming_not_backed_up\2016-01-28_antibiotics\asc810_tetracycline_ON30c_1040\copy-relevant-fluor-images\F\analysis\complete_analysis_16-Feb-2016.mat

myFPPrefix = 'fluor_';
myPCPrefix = '_'; % does not exist
myTresholdPercentile = 97;

myFolder='F:\A_Tans1_step1_incoming_not_backed_up\2016-01-28_antibiotics\asc810_tetracycline_ON30c_1040\copy-relevant-fluor-images\D\';
fluor_checkinglvls_v2

myFolder='F:\A_Tans1_step1_incoming_not_backed_up\2016-01-28_antibiotics\asc810_tetracycline_ON30c_1040\copy-relevant-fluor-images\E\';
fluor_checkinglvls_v2

myFolder='F:\A_Tans1_step1_incoming_not_backed_up\2016-01-28_antibiotics\asc810_tetracycline_ON30c_1040\copy-relevant-fluor-images\F\';
fluor_checkinglvls_v2

myTET = [0.33, 0.16, 0.08];
%}

figure(100); clf; hold on;

for i=1:numel(databaseValuesSignalNoise)
    plot(ones(1,numel(databaseValuesSignalNoise{i}))*myTET(i),...
        databaseValuesSignalNoise{i},'x','Color','r','MarkerSize',15,'LineWidth',2);    
end

plot(myTET,databaseValuesMeanSignalNoise,'o','MarkerSize',15,'LineWidth',2,'Color','k');

%set(gca,'xscale','log');

% Set title etc.
title([SCRIPTNAME 10 'Fluor signal in cells exposed to TET'],'Interpreter','None')
MW_makeplotlookbetter(15);
xlabel('Concentration tetracyclin [mM]');
ylabel('Fluor signal [normalized w. background]');

xlim([0, 0.5])
ylim([0, 3]);


