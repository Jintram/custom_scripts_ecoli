

%% Plotting ribosomal data 

%% 
PLOTOUTDIR = 'U:\ZZ_EXPERIMENTAL_DATA\A_Step5_Data_per_project\Ribosomes\someplots\';
SELECTIONFIELD = 'groupID';
DATASETSTOPLOT = {'prrsaMCerulean_pn25PVenus'};
TOPLOTFIELDNAMES = {{'concentrationCorrData', 'rateCorrData','concentrationDualCrossCorrData', 'rateDualCrossCorrData'},...
                    {'concentrationCorrData', 'rateCorrData'}}; 

%%

for DUALCOLORINDEX=1:2
for plotFieldIndex=1:numel(TOPLOTFIELDNAMES{DUALCOLORINDEX})

    disp('-------------------------------------------');
    disp(['Plotting for DUALCOLORINDEX=' num2str(DUALCOLORINDEX) ', TOPLOTFIELDNAMES=' TOPLOTFIELDNAMES{DUALCOLORINDEX}{plotFieldIndex}]);
    
    % Set plotfield
    TOPLOTFIELDNAME=TOPLOTFIELDNAMES{DUALCOLORINDEX}{plotFieldIndex};
    
    % Run plotting script
    plottingGeneralDynamicData
    
end
end
    
%%

winopen(PLOTOUTDIR);
disp('One script to rule them all done.');