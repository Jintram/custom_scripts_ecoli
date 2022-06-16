
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
ourSettings.mypathname = ['F:\A_Tans1_step1_incoming_not_backed_up\' crossCorrData(CROSSCORRINDEX).dateExperiment];
ourSettings.myID       = crossCorrData(CROSSCORRINDEX).groupID;
ourSettings.myGroupID  = crossCorrData(CROSSCORRINDEX).groupID;

% PARAMETERS TO SAVE?
ourSettings.fitTimeCrosscorr = crossCorrData(CROSSCORRINDEX).params.fitTime;
ourSettings.fitTimeMu        = crossCorrData(CROSSCORRINDEX).params.fitTime;

ourSettings.timeFieldName    = crossCorrData(CROSSCORRINDEX).concentrationFieldNames{1};     
ourSettings.fluorFieldName   = crossCorrData(CROSSCORRINDEX).concentrationFieldNames{2};
ourSettings.muFieldName      = crossCorrData(CROSSCORRINDEX).concentrationFieldNames{3};

ourSettings.timeFieldNameDerivative  = crossCorrData(CROSSCORRINDEX).rateFieldNames{1};
ourSettings.fluorDerivativeFieldName = crossCorrData(CROSSCORRINDEX).rateFieldNames{2};
ourSettings.muFieldNameDerivative    = crossCorrData(CROSSCORRINDEX).rateFieldNames{3};

ourSettings.badSchnitzes = crossCorrData(CROSSCORRINDEX).concentrationBadSchnitzes

ourSettings.alreadyRemovedInMatFile = 0;

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

p.dateDir = [ourSettings.mypathname '\'];

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
    ourSettings.badSchnitzes;

crossCorrData(CROSSCORRINDEX).concentrationautoCorrData = ...
    output.concentrationautoCorrData;
crossCorrData(CROSSCORRINDEX).concentrationautoFieldNames = ...
    output.concentrationautoFieldNames;
crossCorrData(CROSSCORRINDEX).concentrationautoBadSchnitzes = ...
    ourSettings.badSchnitzes;

crossCorrData(CROSSCORRINDEX).rateautoCorrData = ...
    output.rateautoCorrData;
crossCorrData(CROSSCORRINDEX).rateautoFieldNames = ...
    output.rateautoFieldNames;
crossCorrData(CROSSCORRINDEX).rateautoBadSchnitzes = ...
    ourSettings.badSchnitzes;

disp(['Done updated crossCorrData ' num2str(CROSSCORRINDEX) '.']);