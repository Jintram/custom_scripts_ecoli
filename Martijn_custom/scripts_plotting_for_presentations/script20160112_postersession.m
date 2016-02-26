

%% register data
% PROCEED WITH CAUTION, as central project file will be edited!!

% When you have analyzed a dataste, code below can be used to add to
% database.
% PROCEED WITH CAUTION, as central project file will be edited!!
%{
CONFIGFILE = 'U:\ZZ_EXPERIMENTAL_DATA\A_Step5_Data_per_project\CRPcAMP\config_projectCRPcAMP.m'
savingFluorDynamicsData
%}

%% make plot

close all;

CONFIGFILE = 'U:\ZZ_EXPERIMENTAL_DATA\A_Step5_Data_per_project\CRPcAMP\config_projectCRPcAMP.m'
SELECTIONFIELD = 'groupID'

%{
DATASETSTOPLOT = {'ASC852-WtCRP-rCRP' 'ASC853-WtCRP-rS70'}
DATASETSTOPLOT = {'WT_plRCRP-GFP'    'WT_plRs70-GFP'}
DATASETSTOPLOT = {'WT_plRCRP-GFP'    'WT_plRs70-GFP' 'dcAMP-extracell800-CRP' 'dcAMP-extracell800-s70'}    
%}

% WT CRP 
%DATASETSTOPLOT = {'ASC852-WtCRP-rCRP' 'WT_plRCRP-GFP'}
% WT s70
%DATASETSTOPLOT = {'ASC853-WtCRP-rS70' 'WT_plRs70-GFP'}
% cAMP CRP
DATASETSTOPLOT = {'dcAMP-extracell800-CRP'}
% cAMP s70
%DATASETSTOPLOT = {'dcAMP-extracell800-s70'}    
% all relevant CRP
%DATASETSTOPLOT = {'ASC852-WtCRP-rCRP' 'WT_plRCRP-GFP' 'dcAMP-extracell800-CRP'}
% all relevant s70
%DATASETSTOPLOT = {'ASC853-WtCRP-rS70' 'WT_plRs70-GFP' 'dcAMP-extracell800-s70'}
% delta cAMP, all
DATASETSTOPLOT = {'dcAMP-extracell800-CRP', 'dcAMP-extracell800-s70'}

%TOPLOTFIELDNAME = 'concentrationCorrData';
TOPLOTFIELDNAME = 'rateCorrData';
%TOPLOTFIELDNAME = 'growthautoData'
%TOPLOTFIELDNAME = 'concentrationautoCorrData'
%TOPLOTFIELDNAME = 'rateautoCorrData'
    
plottingGeneralDynamicData











