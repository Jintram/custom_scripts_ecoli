

%% Plotting ribosomal data 

%% 
PLOTOUTDIR = 'U:\ZZ_EXPERIMENTAL_DATA\Data_per_project\Ribosomes\someplots\';
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

%% ========================================================================
% PART III: Kiviet style plots
% =========================================================================
%
% The script below is "stolen" from the CRP cross-corr script and these two
% scripts should be integrated.
%
%
%CONFIGFILE='U:\ZZ_EXPERIMENTAL_DATA\A_Step5_Data_per_project\Ribosomes\config_project_Ribosomes.m'

% Setting up
global crossCorrData;

SUBDIR = 'overview\';
SELECTIONFIELD = 'groupID';

if ~exist([PLOTOUTDIR SUBDIR],'dir')
    mkdir([PLOTOUTDIR SUBDIR]);
end

%% Loading
if ~exist('STOPRELOADING','var')
    plottingGeneral_sub_loaddataset
    STOPRELOADING = 1;
end

%% Set parameters
if ~exist('GROUPSTOPLOT','var')
    GROUPSTOPLOT={'prrsaMCerulean_pn25PVenus'};
end
if ~exist('FILENAMESperGroupPartI','var')
    % should match GROUPSTOPLOT
    FILENAMESperGroupPartI = {'rrsa_'};
end
if ~exist('FILENAMESperGroupPartII','var')
    % used to distinguis fluor colors 1, 2, etc
    FILENAMESperGroupPartII = {'rrna-cfp_','pn25-yfp_'}; 
end

if ~exist('TOPLOTFIELDNAMES','var')
    % no need to change
    TOPLOTFIELDNAMES    = {'concentrationCorrData', 'rateCorrData'}; 
end

if ~exist('TITLESforgroupsPartI','var')    
    TITLESforgroupsPartI = {''};
end
if ~exist('TITLESforgroupsPartII','var')
    TITLESforgroupsPartII = {'rRNA reporter','pn25 reporter'};
end

Line1Colors = linspecer(numel(GROUPSTOPLOT));
Line2Colors = ones(numel(GROUPSTOPLOT),3)*0; % also allows gray

%% Do it using loop

optionsStruct=struct;

for groupIdx = 1:numel(GROUPSTOPLOT)
    for dualColorIdx = 1:2
        
        % parameters that are used by plotting script
        DATASETSTOPLOT = {GROUPSTOPLOT{groupIdx} GROUPSTOPLOT{groupIdx}};
        DUALCOLORINDICES =  [dualColorIdx dualColorIdx];
        LINECOLORS = {Line1Colors(groupIdx,:), Line2Colors(groupIdx,:)};

        % plotting script
        optionsStruct.STOPRELOADING=1;
        [hCC,output]=plottingGeneral_v2_CCs(DATASETSTOPLOT,DUALCOLORINDICES,LINECOLORS,SELECTIONFIELD,TOPLOTFIELDNAMES,optionsStruct)        

        % plot cosmetics
        figure(hCC.Number); 
        title([TITLESforgroupsPartI{groupIdx},TITLESforgroupsPartII{dualColorIdx}]);
        xlim([-10*output.hrsPerDoublingMean,10*output.hrsPerDoublingMean]);
        
        % Save file adjusted xlims
        for extIdx = 1:numel(EXTENSIONS)
            fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_xlimsadj_' FILENAMESperGroupPartI{groupIdx} FILENAMESperGroupPartII{dualColorIdx} '.' EXTENSIONS{extIdx}];
            if strcmp(EXTENSIONS{extIdx},'eps'), saveas(hCC,fileName,'epsc'); else saveas(hCC,fileName); end
        end
        
        xlim([-10,10]);
        
        % Save file with same xlims
        for extIdx = 1:numel(EXTENSIONS)
            fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_' FILENAMESperGroupPartI{groupIdx} FILENAMESperGroupPartII{dualColorIdx} '.' EXTENSIONS{extIdx}];
            if strcmp(EXTENSIONS{extIdx},'eps'), saveas(hCC,fileName,'epsc'); else saveas(hCC,fileName); end
        end
    end
end

disp('Done plotting Cross corrs');
winopen([PLOTOUTDIR SUBDIR]);