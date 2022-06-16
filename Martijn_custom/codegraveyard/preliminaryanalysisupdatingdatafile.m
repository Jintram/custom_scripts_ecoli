%% It is very important to supply the path to the correct dataset here!
CONFIGFILEPATH = 'F:\A_Tans1_step1_incoming_not_backed_up\2015-10-20\Schnitzcells_Analysis_Config_2015_10_20_pos8.xlsx';
ANALYSISTYPE = 1; % 1 = preliminary, 2 = full

% Parameters you SHOULD NOT change
EXCELREADSTART = 14; % line where list of parameters starts in Excel file.

% Read configuration settings from excel file.

% One can execute this section again to reload settings 
% Note that these are then not immediate parsed to the "p" struct.
[ourSettings, alldata] = MW_readsettingsfromexcelfile(CONFIGFILEPATH)

% Double check whether to continue
% This step is mainly to prevent accidentally running whole script 
analysisTypes = {'preliminary','full'}

myAnswer = questdlg(['Loaded ' CONFIGFILEPATH ', for ' analysisTypes{ANALYSISTYPE} ' analysis do you want to continue?'],'Confirmation required.','Yes','No','No');
if strcmp(myAnswer, 'No') || strcmp(myAnswer, 'Cancel')
    error('Analysis aborted.');
end

% Now make a vector with parameter values that schnitzcells scripts can handle. (and some misc. other admin)

disp('Now creating ''p'' struct from ourSettings struct.');

% Create the p vector which holds all parameters and is fed into, and also
% returned by most functions of the schnitzcells analysis software.
% Use "ourSettings" struct as a base.
% DJK_initschnitz only checks for errors and adds a few parameters based on
% the already given paramteres.
% ===
% TODO this can be done more elegantly (but note that existence of two
% vectors, "ourSettings" and "p", allows user to update "ourSettings" vector 
% intermediately.
p = DJK_initschnitz(ourSettings.positionName,ourSettings.movieDate,'e.coli.amolf','rootDir',...
    ourSettings.rootDir, 'cropLeftTop',ourSettings.cropLeftTop, 'cropRightBottom',ourSettings.cropRightBottom,...
    'fluor1',ourSettings.fluor1,'fluor2',ourSettings.fluor2,'fluor3',ourSettings.fluor3,...
    'setup',ourSettings.setup,'softwarePackage',ourSettings.softwarePackage,'camera',ourSettings.camera)

% Manually make sure image dir is correct
% (This is done to accomodate cropping.)
p.imageDir = [ourSettings.rootDir ourSettings.movieDate '\' ourSettings.positionName '\']

% Set framerange according to analysis type
if ANALYSISTYPE == 1 % fast analysis
    currentFrameRange = ourSettings.frameRangePreliminary;
elseif ANALYSISTYPE == 2 % full analysis
    currentFrameRange = ourSettings.frameRangeFull;
else
    error('No analysis type speficied');
end


% =========================================================================
% (After loading ourSettings from excel file.)
p = DJK_initschnitz([ourSettings.positionName ourSettings.cropSuffix],ourSettings.movieDate,'e.coli.amolf','rootDir',...
    ourSettings.rootDir, 'cropLeftTop',ourSettings.cropLeftTop, 'cropRightBottom',ourSettings.cropRightBottom,...
    'fluor1',ourSettings.fluor1,'fluor2',ourSettings.fluor2,'fluor3',ourSettings.fluor3,...
    'setup',ourSettings.setup,'softwarePackage',ourSettings.softwarePackage,'camera',ourSettings.camera)
% =========================================================================



    % Load this dataset using
    % =========================================================================
    p = DJK_initschnitz([ourSettings.positionName ourSettings.cropSuffix], ourSettings.movieDate,'e.coli.AMOLF','rootDir',...
     ourSettings.rootDir, 'cropLeftTop', ourSettings.cropLeftTop, 'cropRightBottom', ourSettings.cropRightBottom,...
         'fluor1',ourSettings.fluor1,...
         'fluor2',ourSettings.fluor2,...
         'fluor3',ourSettings.fluor3,...
         'setup',ourSettings.setup,...
         'softwarePackage',ourSettings.softwarePackage,...
         'camera',ourSettings.camera);

    [p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);
    % =========================================================================
    
    
    
if ANALYSISTYPE==1 % fast

    % Output directory
    myOutputDir = [p.dateDir 'outputSummary\'];

    % Obtain schnitzcells
    [p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

    % Fit mu
    [fitTime, fitMu] = DJK_analyzeMu(p, schnitzcells, 'onScreen', 1,'fitTime',[0 10000],'DJK_saveDir',myOutputDir);
    close(gcf);
    %[fitTime, fitMu] = DJK_analyzeMu(p, schnitzcells, 'xlim', [0 900], 'onScreen', 1,'fitTime',[0 10000]);

    % Again loop over colors, plot signal for each
    fluorColors = {'fluor1','fluor2','fluor3'}; % ugly but compatible w. legacy - MW
    fitFluorMean = nan(1,3); fitFluorVariance = nan(1,3);
    for colorIdx = 1:3

        % Adminstration
        % Select current string w. color identifier (fluor1, ..)    
        currentFluor = fluorColors{colorIdx};
        if strcmp(p.(currentFluor),'none') || strcmp(p.(currentFluor),'None')
            disp([currentFluor ' not set, assuming you don''t have more fluor colors.']);
            break
        end

        % plot fluor behavior
        fluorFieldName = [upper(p.(currentFluor)(1)) '5_mean_all'];
        [fitFluorMean(colorIdx), fitFluorVariance(colorIdx)] = DJK_plot_avColonyOverTime(p, schnitzcells, fluorFieldName, 'fitTime', fitTime, 'onScreen', 1,'DJK_saveDir',myOutputDir);
        close(gcf);
        %DJK_plot_avColonyOverTime(p, schnitzcells, 'C5_mean_all', 'xlim', [0 900], 'ylim', [0 150], 'fitTime', fitTime, 'onScreen',1);
    end

    numberofSchnitzesInLast = MW_countcellsinlastframe(schnitzcells);

    % Write to summary file
    summaryParametersNames = {'#Schnitzes in last frame', 'fitMu', 'fitTime(1)', 'fitTime(2)', 'fitFluorMean1', 'fitFluorVariance1','fitFluorMean2', 'fitFluorVariance2','fitFluorMean3', 'fitFluorVariance3'}
    summaryParameters = [numberofSchnitzesInLast, fitMu, fitTime, fitFluorMean, fitFluorVariance]

    % obtain the number of this position (e.g. pos1 => 1)
    posNumber = str2num(p.movieName(regexp(p.movieName,'[*\d]')))
    
    lineToWriteTo = num2str(posNumber+1);
    
    % Output to excel sheet
    xlswrite([myOutputDir 'summaryParametersPreliminary.xls'],{'Identifier'},['A1:A1'])
    xlswrite([myOutputDir 'summaryParametersPreliminary.xls'],{[p.movieDate '_' p.movieName]},['A' lineToWriteTo ':' 'A' lineToWriteTo])
    
    for i = 1:numel(summaryParameters)
        letterToWriteTo = ['' i+64+1]; % +1 since start at B
        %disp(['writing ' [letterToWriteTo lineToWriteTo]]);
        xlswrite([myOutputDir 'summaryParametersPreliminary.xls'],{summaryParametersNames{i}},[letterToWriteTo '1:' letterToWriteTo '1'])
        xlswrite([myOutputDir 'summaryParametersPreliminary.xls'],summaryParameters(i),[letterToWriteTo lineToWriteTo ':' letterToWriteTo lineToWriteTo])
    end

    % Save output to .mat file (update the matfile)
    if exist([myOutputDir 'summaryParametersPreliminary.mat'],'file') == 2
        load([myOutputDir 'summaryParametersPreliminary.mat'],'thedata');
    end
    thedata(posNumber).summaryParameters = summaryParameters;
    thedata(posNumber).ourSettings = ourSettings;
    thedata(posNumber).p = p; % "ourSettings" and "p" are a bit redundant for historic reasons.
    save([myOutputDir 'summaryParametersPreliminary.mat'],'thedata','summaryParametersNames');
    
disp('Done making summary preliminary analysis.');    
    
end