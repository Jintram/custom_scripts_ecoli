
% This script was made to re-run the main analysis for data sets taken
% earlier, that were analyzed at a time the analysis was different and less
% extensive.
%
% General approach to run this file:
% - Load crossCorrData.mat (use plottingGeneralDynamicData)
% - Set CROSSCORRINDEX and MYFLUOR appropriately
% - Run script
% - Save crossCorrData.mat using appropriate lines from savingFluorDynamicsData
%
% MW, 2015/02/29

% PARAMETERS to set
CROSSCORRINDEX=17
MYFLUOR = {'g','',''}; % fluors used, cells correspond to p.fluor1, etc

% for settings
% ===
settings.mypathname = ['F:\A_Tans1_step1_incoming_not_backed_up\' crossCorrData(CROSSCORRINDEX).dateExperiment];
settings.myID       = crossCorrData(CROSSCORRINDEX).groupID;
settings.myGroupID  = crossCorrData(CROSSCORRINDEX).groupID;

% PARAMETERS TO SAVE?
settings.fitTimeCrosscorr = crossCorrData(CROSSCORRINDEX).params.fitTime;
settings.fitTimeMu        = crossCorrData(CROSSCORRINDEX).params.fitTime;

settings.timeFieldName    = crossCorrData(CROSSCORRINDEX).concentrationFieldNames{1};     
settings.fluorFieldName   = crossCorrData(CROSSCORRINDEX).concentrationFieldNames{2};
settings.muFieldName      = crossCorrData(CROSSCORRINDEX).concentrationFieldNames{3};

settings.timeFieldNameDerivative  = crossCorrData(CROSSCORRINDEX).rateFieldNames{1};
settings.fluorDerivativeFieldName = crossCorrData(CROSSCORRINDEX).rateFieldNames{2};
settings.muFieldNameDerivative    = crossCorrData(CROSSCORRINDEX).rateFieldNames{3};

settings.badSchnitzes = crossCorrData(CROSSCORRINDEX).concentrationBadSchnitzes

settings.alreadyRemovedInMatFile = 0;

% for p
% ===
p.movieDate  = crossCorrData(CROSSCORRINDEX).dateExperiment;

% extract name
identifierString = crossCorrData(CROSSCORRINDEX).identifier;
underscoreidx = strfind(identifierString,'_');
movieName = identifierString(underscoreidx+1:end);

p.movieName  = movieName;

p.fluor1 = MYFLUOR{1};
p.fluor2 = MYFLUOR{2};
p.fluor3 = MYFLUOR{3};

p.dateDir = [settings.mypathname '\'];

% now rerun the analysis
% ===
runsections = 'rerunfullanalysis';
MW_analysis_attempt2_matlabinsteadexcel_plusGUI
    


%% now update crossCorrData

crossCorrData(CROSSCORRINDEX).growthautoData = ...
    output.growthautoCorrData;
crossCorrData(CROSSCORRINDEX).growthautoFieldNames = ...
    output.growthautoFieldNames;
crossCorrData(CROSSCORRINDEX).growthautoBadSchnitzes = ...
    settings.badSchnitzes;

crossCorrData(CROSSCORRINDEX).concentrationautoCorrData = ...
    output.concentrationautoCorrData;
crossCorrData(CROSSCORRINDEX).concentrationautoFieldNames = ...
    output.concentrationautoFieldNames;
crossCorrData(CROSSCORRINDEX).concentrationautoBadSchnitzes = ...
    settings.badSchnitzes;

crossCorrData(CROSSCORRINDEX).rateautoCorrData = ...
    output.rateautoCorrData;
crossCorrData(CROSSCORRINDEX).rateautoFieldNames = ...
    output.rateautoFieldNames;
crossCorrData(CROSSCORRINDEX).rateautoBadSchnitzes = ...
    settings.badSchnitzes;

disp(['Done updated crossCorrData ' num2str(CROSSCORRINDEX) '.']);