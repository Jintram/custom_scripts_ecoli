


%% A collective script that makes plots for the miscelaneous parts of the thesis.

OUTPUTFOLDER='\\storage01\data\AMOLF\users\wehrens\Latex3\Thesis\Chapter2_Methods\Figures\MatlabExport\';

%% Some plots relating to the weighing of cross correlations for technical chapter

general_weighinCCs_setup

figure(hInput)
% Save it
extensions={'svg','tif','fig'};
for extensionIdx=1:3, extension=extensions{extensionIdx};
    saveas(hInput,[OUTPUTFOLDER extension '_' 'Technical_CC_weighing_A' '_.' extension]);
end

figure(hCC)
% Save it
extensions={'svg','tif','fig'};
for extensionIdx=1:3, extension=extensions{extensionIdx};
    saveas(hCC,[OUTPUTFOLDER extension '_' 'Technical_CC_weighing_B' '_.' extension]);
end

