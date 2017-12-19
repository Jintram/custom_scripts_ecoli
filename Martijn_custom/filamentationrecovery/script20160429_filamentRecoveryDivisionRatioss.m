
%% 

% NOTE ABOUT CorrectedQuickDiv
% ===
% Note that I tried to correct for cells that divide fast, but that this is
% impossible, since we don't know the orientation of the cells.

% Plotting division ratios
if ~exist('WHATDATA','var')
    WHATDATA = 'sulA';
    WHATDATA = 'temperature';
    %WHATDATA = 'simulated';
    WHATDATA = 'tetracycline';
end

WINDOWBORDERS   = [3:6:36];
%WINDOWBORDERS   = [1.8:5.2:28];
BARCOLORS       = [240,155,34; 45 177 65; 37 156 190; 131 84 162; 241 88 58]./255;

USESYMMETRY=1;
HISTNRBINS=50;

FIGURENUMBERS=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % TODO this option does not work any more

%PLOTCOLORS = [0 0 1;0 0 1;0 0 1];
%PLOTCOLORS = [0,0,139 ; 0,191,255  ; 	135,206,250]./255; % shades of blue
%PLOTCOLORS = [255, 0, 0; 204, 0, 204 ; 255, 102, 0]./255; % shades of blue
%Later: PLOTCOLORS = linspecer(numel(datasetsPaths));

if ~exist('LENGTHFIELD','var')
    %LENGTHFIELD = 'areaPixels';
    LENGTHFIELD = 'length_fitNew';
    %LENGTHFIELD = 'length_skeleton';
    %LENGTHFIELD = 'cellLengths';
end

TIMEFIELD = 'time';
%TIMEFIELD = 'times';

if ~exist('PLOTSAVEDIR','var') & ~exist('NOSAVEPLEASE','var');
    error('Set PLOTSAVEDIR please. Or set NOSAVEPLEASE');
    % PLOTSAVEDIR = 'D:\Local_Data\Dropbox\Dropbox\Filamentation recovery\MW\figures_new\matlab_export\';
end
if ~exist('NOSAVEPLEASE','var')
    NOSAVEPLEASE=0;
end

if ~exist('SPECIALCASE','var')
    SPECIALCASE=0;
end

% user given parameters    
LEFTX =  0;   

%{
% Please execute this code to define plot dir and make log of scriptname
PLOTSAVEDIR = 'U:\PROJECTS\B_filamentationRecovery\Contributions_Martijn\figures_fig_files\';
PLOTSAVEDIR='D:\Local_Data\Dropbox\Dropbox\Filamentation recovery\MW\figures_new\matlab_export\'
makenoteinreadme(PLOTSAVEDIR,'script20160429_filamentRecoveryDivisionRatioss');
%}

%% 
if ~exist('RUNSECTIONSFILADIV','var')
    disp('Set RUNSECTIONSFILADIV to run this script, set RUNSECTIONSFILADIV to ''all'' to run everything.');
end

%% The dataset paths + specific information

if strcmp(WHATDATA, 'sulA')
    datasetsPaths = ...
        { ...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-08_FilaRecovery_sulA_recovery_200uM_IPTG\pos1crop\data\pos1crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-08_FilaRecovery_sulA_recovery_200uM_IPTG\pos2crop\data\pos2crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-08_FilaRecovery_sulA_recovery_200uM_IPTG\pos3crop\data\pos3crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-08_FilaRecovery_sulA_recovery_200uM_IPTG\pos4crop\data\pos4crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-08_FilaRecovery_sulA_recovery_200uM_IPTG\pos7crop\data\pos7crop-Schnitz.mat'...
        }
    switchTimes = [0 0 0 0 0];
    SPECIALCASE=0;
elseif strcmp(WHATDATA, 'temperature')
    datasetsPaths = ...
        { ...
        'G:\EXPERIMENTAL_DATA_2016\2016-03-23_FilaRecovery_asc777_42C\pos4crop\data\pos4crop-Schnitz.mat',...
            ... No fluor signal
        ['G:\EXPERIMENTAL_DATA_2016\2016-04-07_FilaRecovery_asc777\pos2crop\data\pos2crop-Schnitz.mat'],...
            ... 
        }
    switchTimes = [450, 329];
    SPECIALCASE=0;
elseif strcmp(WHATDATA, 'deltaMinTET')
    datasetsPaths = ...
        { ...
        'H:\EXPERIMENTAL_DATA_2017\2017-09-22_FilaRecovery_asc1035_DeltaMinCDE\pos1cropb\data\pos1cropb-Schnitz.mat',...        
        'H:\EXPERIMENTAL_DATA_2017\2017-09-22_FilaRecovery_asc1035_DeltaMinCDE\pos2cropa\data\pos2cropa-Schnitz.mat',...
        'H:\EXPERIMENTAL_DATA_2017\2017-09-22_FilaRecovery_asc1035_DeltaMinCDE\pos2cropb\data\pos2cropb-Schnitz.mat',...
        'H:\EXPERIMENTAL_DATA_2017\2017-09-22_FilaRecovery_asc1035_DeltaMinCDE\pos2cropc\data\pos2cropc-Schnitz.mat',...
        'H:\EXPERIMENTAL_DATA_2017\2017-09-22_FilaRecovery_asc1035_DeltaMinCDE\pos3cropa\data\pos3cropa-Schnitz.mat',...
            ...         
            ... 
        }
    switchTimes = [5 5 5 5 5];
        % frame 600=4:26, switch @4:26 (but 5min delay)
    SPECIALCASE=0;
elseif strcmp(WHATDATA, 'tetracycline')
    
    % Load data from gathered file:
    %SPECIALCASE=1;
    %load('G:\FilamentationRecoveryData\Dmitry\3_Rutgers_project\Recent Data\rutgers old data\DE analysis\analysis of clusters\all_rutgers_divisions_data_gen_shifted.mat');    
    %datasetsPaths=exp_names_all;
    
    % Or from Rutger's files
    ONESTOTAKE=[1:5]; % ONESTOTAKE=[1:5]; ONESTOTAKE=[6:8]; ONESTOTAKE=[9:11];
    if exist('SPECIALONESTOTAKE','var')
        ONESTOTAKE=SPECIALONESTOTAKE;
        %clear SPECIALONESTOTAKE;
    end
    disp(['Note that certain datasets are selected (' num2str(ONESTOTAKE) ').']);
    switchTimes = ...
        ...[1373.5417      587.66667      261.13333        741.325      252.38333      987.86667        1017.95         997.05      881.90833      1018.7083         1394.6]; % based on 50% max growth rate
        ... [1.2156    0.5064    0.0512    0.6032    0.0948    0.8585    0.7400    0.7954    0.7673    0.8160    1.2737]*1.0e+03; % based on 1/10th max. growth rate
        [890.9800  404.7500         0  529.7600  0 ... % based on calculated switch times
            0 0 0 0 0 0]; % NOTE THAT HERE NO TIMES ARE DETERMINED!!
        %warning('change this back!!');
        % switchTimes determined by 
        % script20161222_branchplottingrutgerdata
    switchTimes = switchTimes(ONESTOTAKE);
    datasetsPaths = ...
        { ...
        ... Note that the parameter ONESTOANALYZE makes a subselection of this data.
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos3_long.mat',...             1
        ... ^ Comes from ..\2013-12-09\pos3crop\data\ (note switchtimes are in descriptions)
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos4.mat',...                  2
        ... ^ Comes from ..\2013-09-24\pos4crop\data\ (ugly images)
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos4_long.mat',...             3
        ... ^ Comes from ..\2013-12-16\pos4crop\data\ (ugly images)
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos5.mat',...                  4
        ... ^ Comes from ..\2013-09-24\pos5crop\data\ (very little data)
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos5_long.mat',...             5
        ... ^ Comes from ..\2013-12-16\pos5crop\data\ (ugly images)
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\2uM_pos2.mat',...          6
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\2uM_pos4.mat',...          7
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\2uM_pos6.mat',...          8
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\10uM_pos1.mat',...         9
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\10uM_pos3.mat',...         10
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\10uM_pos6_long.mat',...    11  
        }    
    datasetsPaths={datasetsPaths{ONESTOTAKE}};    
elseif strcmp(WHATDATA, 'tetracycline_skel')
    % This option of WHATDATA is just for double checking whether indeed 
    
    switchTimes = ...
        ...[1373.5417      587.66667      261.13333        741.325      252.38333      987.86667        1017.95         997.05      881.90833      1018.7083         1394.6]; % based on 50% max growth rate
        ... [1.2156    0.5064    0.0512    0.6032    0.0948    0.8585    0.7400    0.7954    0.7673    0.8160    1.2737]*1.0e+03; % based on 1/10th max. growth rate
        [890.9800  404.7500         0  529.7600  0]; % based on calculated switch times
        %warning('change this back!!');
        % switchTimes determined by 
        % script20161222_branchplottingrutgerdata
    
    datasetsPaths = ...
        { ...
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-09\pos3crop\data\pos3crop-Schnitz.mat',...
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos4crop\data\pos4crop-Schnitz.mat',...
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos4crop\data\pos4crop-Schnitz.mat',...
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos5crop\data\pos5crop-Schnitz.mat',...
        ...'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos5crop\data\pos5crop-Schnitz.mat',...
        };
elseif strcmp(WHATDATA, 'simulated')
    % now there already should be a simulated schnitzcells, use that one ..
    if exist('simulatedschnitzcells','var')
        schnitzcells=simulatedschnitzcells;
    % or use old stored one
    else
        datasetsPaths = ...
            { ...
            'D:\Local_Software\Martijn_extensions\Martijn_custom\filamentationrecovery\simulatedSchnitzcells\schnitz1_naivemodel_fastrecharge.mat', ...
            'D:\Local_Software\Martijn_extensions\Martijn_custom\filamentationrecovery\simulatedSchnitzcells\schnitz1_naivemodel_longrecharge.mat' ...
            };
    end
else
    error('No data loaded..');
end

PLOTCOLORS = linspecer(numel(datasetsPaths));

%% Gather data
% Loop over datasets, then over individual schnitzes. Create a parameter
% structure that looks as follows:
%
% - myLengthNewborns{dataSetIndex}(schnitz) 
% - myLengthParents{dataSetIndex}(schnitz)
% - myLengthSumNewborns{dataSetIndex}(schnitz)
%
% Where myLengthNewborns gives the length of a newborn, and myLengthParents
% gives the length of the corresponding parent. Because there are some
% subtleties (growth, geometrical distortion poles), myLengthSumNewborns is
% a better measure of the total length before division, and is simply the
% summed length of the two daughter cells.
if any(strcmp(RUNSECTIONSFILADIV,{'all','loadData'}))

    disp('Loading data..');
    schnitzcount=0;

    myLifeTimeParents = {}; 
    myLifeTimesSchnitzes = {};
    myNewBornSchnitzNrs = {};
    allLengths = {};
    birthTimes = {};
    myLengthNewborns = {};
    myLengthSumNewborns = {};
    %myLengthParentsCorrectedQuickDiv={};
    %myLifeTimeParentsCorrectedQuickDiv={};
    %myLengthSumNewbornsCorrectedQuickDiv={};
    %listWhichCorrectedQuickDiv={};
    myLengthSumNewborns = {};
    myCellLengthsPerSchnitz = {}; 
    myAddedLengthPerSchnitz = {};
    usedSchnitzes = {};

    for dataSetIndex = 1:numel(datasetsPaths)

        if ~exist('simulatedschnitzcells','var')
            if ~SPECIALCASE
                schnitzcells = loadandrename(datasetsPaths{dataSetIndex});
            else
                schnitzcells=S_all_shifted{dataSetIndex};
            end
        end
        schnitzcount=schnitzcount+numel(schnitzcells);

        %% 
        figureIndex=FIGURENUMBERS(dataSetIndex);

        %%


        % Finding parent with each daughter
        %{
        myLengthNewborns = []; myLengthParents = [];
        for i = 1:numel(schnitzcells)

            LengthNewborn = schnitzcells(i).(LENGTHFIELD)(1);

            parentSchnitz = schnitzcells(i).P;    

            if parentSchnitz ~=0

                LengthParent = schnitzcells(parentSchnitz).(LENGTHFIELD)(end);

                myLengthNewborns(end+1) =   LengthNewborn;
                myLengthParents(end+1) =    LengthParent;

            end

        end
        %}

        %% Finding 2 daughters with each parent
        myLengthNewborns{dataSetIndex} = []; myLengthParents{dataSetIndex} = []; myLengthSumNewborns{dataSetIndex}=[];
        myNewBornSchnitzNrs{dataSetIndex} = []; %myDaughter1schnitzNrs{dataSetIndex} = []; myDaughter2schnitzNrs{dataSetIndex} = [];
        myLifeTimeParents{dataSetIndex}=[];
        usedSchnitzes{dataSetIndex} = zeros(1,numel(schnitzcells));
        %myLengthParentsCorrectedQuickDiv{dataSetIndex}=[];
        %myLifeTimeParentsCorrectedQuickDiv{dataSetIndex}=[];
        %myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}=[];
        %listWhichCorrectedQuickDiv{dataSetIndex}=[];
        
        for binIdx = 1:numel(schnitzcells) % todo: rename iteration parameter; incorrect name --> schnitzIdx

            LengthParent = schnitzcells(binIdx).(LENGTHFIELD)(end);

            daughterSchnitz1 = schnitzcells(binIdx).D;
            daughterSchnitz2 = schnitzcells(binIdx).E;                                    
            
            if ~any([daughterSchnitz1,daughterSchnitz2]==0) & ...
                    ~isempty(schnitzcells(daughterSchnitz1).(LENGTHFIELD)) & ...
                    ~isempty(schnitzcells(daughterSchnitz2).(LENGTHFIELD)) 
                
                % create administration which schnitzes were used for this
                % calculation (last/dead ends are pro'lly excluded)
                usedSchnitzes{dataSetIndex}(daughterSchnitz1) = 1;
                usedSchnitzes{dataSetIndex}(daughterSchnitz2) = 1;
                
                % Then calculate ratios
                lengthDaughterSchnitz1 = schnitzcells(daughterSchnitz1).(LENGTHFIELD)(1);
                lengthDaughterSchnitz2 = schnitzcells(daughterSchnitz2).(LENGTHFIELD)(1);            

                % daughter 1
                myLengthNewborns{dataSetIndex}(end+1) =     lengthDaughterSchnitz1;
                myLengthParents{dataSetIndex}(end+1) =      LengthParent;
                myLengthSumNewborns{dataSetIndex}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2;
                myNewBornSchnitzNrs{dataSetIndex}(end+1) =  daughterSchnitz1;
                % daughter 2
                myLengthNewborns{dataSetIndex}(end+1) =     lengthDaughterSchnitz2;
                myLengthParents{dataSetIndex}(end+1) =      LengthParent;
                myLengthSumNewborns{dataSetIndex}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2;
                myNewBornSchnitzNrs{dataSetIndex}(end+1) =  daughterSchnitz2;

                % also note down lifetime of parent
                if isfield(schnitzcells,'interDivTime')
                    myLifeTimeParents{dataSetIndex}(end+1) = schnitzcells(binIdx).interDivTime; % daughter 1
                    myLifeTimeParents{dataSetIndex}(end+1) = schnitzcells(binIdx).interDivTime; % daughter 2
                else
                    myLifeTimeParents{dataSetIndex}(end+1) = NaN; % daughter 1
                    myLifeTimeParents{dataSetIndex}(end+1) = NaN; % daughter 2
                end                               
                
                %{
                % correct the data for young mothers, i.e. the pattern breaks
                % down when cells divide quickly, either because 1 cell can be
                % considered to have 3 daughters, or because the ring cannot
                % rearrange fast enough
                if schnitzcells(i).interDivTime<20

                    % first order correction
                    correctedParentSchnitz = schnitzcells(i).P;
                    % can't correct though if unknown mother
                    if correctedParentSchnitz==0, 
                        correctedParentSchnitz=i; 
                        currentParentLifeTime = schnitzcells(i).interDivTime;
                        lengthCousin=0;                    
                        warning('First cell not corrected'); 
                        listWhichCorrectedQuickDiv{dataSetIndex}(end+1) = 0; listWhichCorrectedQuickDiv{dataSetIndex}(end+1) = 0;
                    else
                        % start correction
                        listWhichCorrectedQuickDiv{dataSetIndex}(end+1) = 1; listWhichCorrectedQuickDiv{dataSetIndex}(end+1) = 1;
                        currentParentLifeTime = schnitzcells(i).interDivTime + schnitzcells(correctedParentSchnitz).interDivTime;                    

                        % find length cousin    
                        potentialCousins = [schnitzcells(correctedParentSchnitz).D schnitzcells(correctedParentSchnitz).E];
                        theCousin = potentialCousins(potentialCousins~=i); % not the mother itself                
                        theLengths = schnitzcells(theCousin).(LENGTHFIELD);
                        lengthCousin = theLengths(end);
                    end

                    % Use this data
                    LengthGrandParent = schnitzcells(correctedParentSchnitz).(LENGTHFIELD)(end);
                    % Store it
                    myLengthParentsCorrectedQuickDiv{dataSetIndex}(end+1) = LengthGrandParent;
                    myLengthParentsCorrectedQuickDiv{dataSetIndex}(end+1) = LengthGrandParent;
                    myLifeTimeParentsCorrectedQuickDiv{dataSetIndex}(end+1) = currentParentLifeTime;
                    myLifeTimeParentsCorrectedQuickDiv{dataSetIndex}(end+1) = currentParentLifeTime;

                    % calculate length with it
                    myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2+lengthCousin;
                    myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2+lengthCousin;
                else
                    % no correction
                    listWhichCorrectedQuickDiv{dataSetIndex}(end+1) = 0; listWhichCorrectedQuickDiv{dataSetIndex}(end+1) = 0;
                    myLengthParentsCorrectedQuickDiv{dataSetIndex}(end+1) = LengthParent;
                    myLengthParentsCorrectedQuickDiv{dataSetIndex}(end+1) = LengthParent;
                    myLifeTimeParentsCorrectedQuickDiv{dataSetIndex}(end+1) = schnitzcells(i).interDivTime;
                    myLifeTimeParentsCorrectedQuickDiv{dataSetIndex}(end+1) = schnitzcells(i).interDivTime;
                    myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2;
                    myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2;
                end
                %}

                % For later reference, make lookup table 
                %myDaughter1schnitzNrs{dataSetIndex}(end+1) = daughterSchnitz1;
                %myDaughter2schnitzNrs{dataSetIndex}(end+1) = daughterSchnitz2;                        
            
            else
                if ~any([daughterSchnitz1,daughterSchnitz2]==0)
                    if isempty(schnitzcells(daughterSchnitz1).(LENGTHFIELD)) | ...
                            isempty(schnitzcells(daughterSchnitz2).(LENGTHFIELD))

                        warning('Skipping schnitz because length not available..');
                    end
                end
            end
            

        end

        %% Gather data on interdivision time vs. lifetime of cell
        birthTimes{dataSetIndex} = NaN(1,numel(schnitzcells));
        allLengths{dataSetIndex} = NaN(1,numel(schnitzcells));
        myCellLengthsPerSchnitz{dataSetIndex} = NaN(1,numel(schnitzcells));
        myAddedLengthPerSchnitz{dataSetIndex} = NaN(1,numel(schnitzcells));
        
        if isfield(schnitzcells,'interDivTime')
            myLifeTimesSchnitzes{dataSetIndex} = [schnitzcells.interDivTime];
        end
        
        for binIdx = 1:numel(schnitzcells)
            
            thisSchnitzLengths = schnitzcells(binIdx).(LENGTHFIELD);
            allLengths{dataSetIndex}(binIdx) = thisSchnitzLengths(1);
            birthTimes{dataSetIndex}(binIdx) = schnitzcells(binIdx).(TIMEFIELD)(1);
            
            % Also additionally create length vs. added size plot
            % (relating to "adder" discussion)
            if ~any([schnitzcells(binIdx).D,schnitzcells(binIdx).E]==0)
                myCellLengthsPerSchnitz{dataSetIndex}(binIdx) = schnitzcells(binIdx).(LENGTHFIELD)(1); 
                myAddedLengthPerSchnitz{dataSetIndex}(binIdx) = schnitzcells(binIdx).(LENGTHFIELD)(end)-schnitzcells(binIdx).(LENGTHFIELD)(1); 
            end

        end

        if ~isfield(schnitzcells,'interDivTime')
            warning('interDivTime field was not found')
        end

    end

    disp('Loading data done.');
    disp(['N=' num2str(schnitzcount) ' cells were loaded.']);
    
end

%% Calculate ratios
if any(strcmp(RUNSECTIONSFILADIV,{'all','calculateRatios'}))

    Ratios={}; 
    for dataSetIndex = 1:numel(datasetsPaths)

        % calculate ratios
        Ratios{dataSetIndex} = myLengthNewborns{dataSetIndex}./myLengthSumNewborns{dataSetIndex};
        %RatiosCorrectedQuickDiv{dataSetIndex} = myLengthNewborns{dataSetIndex}./myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex};
        
    end
    
end

%% Now calculate ratios and plot them
if any(strcmp(RUNSECTIONSFILADIV,{'all','rawplot'}))

    if ~exist('SANITYCHECKSYMMETRY','var')
        SANITYCHECKSYMMETRY=1;
    end

    % RatiosCorrectedQuickDiv = {};
    for dataSetIndex = 1:numel(datasetsPaths)
        
        %% Gathering cosmetic parameters        

        % calculated parameters
        if ~exist('RIGHTX','var')
            if strcmp(LENGTHFIELD,'length_skeleton')
                RIGHTX = 30;        
                LONGESTNORMALDIVSIZEPARENT=5;
            elseif strcmp(LENGTHFIELD,'areaPixels')
                RIGHTX = 15000;
                LONGESTNORMALDIVSIZEPARENT = 2000;
            else
                RIGHTX  = max(myLengthParents{dataSetIndex})*1.1;
                LONGESTNORMALDIVSIZEPARENT=5;
            end
        end

        % create figure
        figure(1); 

        % target length line
        x = (LONGESTNORMALDIVSIZEPARENT):RIGHTX;
        y = (LONGESTNORMALDIVSIZEPARENT/2)./x;                
        
        %% actual plotting
        
        % plot helping lines at 1/2n
        % for i=1:5
        %     plot([rightx, LEFTX], [.5/i .5/i],':','Color',[.5 .5 .5],'LineWidth',2)
        %     plot([rightx, LEFTX], 1-[.5/i .5/i],':','Color',[.5 .5 .5],'LineWidth',2)
        % end
        N=5;
        for binIdx=1:N
            for j = 1:(binIdx*2-1)
                plot([0, RIGHTX], [(j)/(2*binIdx) (j)/(2*binIdx)],'-','Color',[.5 .5 .5],'LineWidth',N-binIdx+1)
            end
        end
        
        % reference child
        plot(x,y,'--','Color','k','LineWidth',2);
        
        % plot length vs. ratios
        plot(myLengthSumNewborns{dataSetIndex},Ratios{dataSetIndex},'o', 'Color', PLOTCOLORS(dataSetIndex,:),'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',PLOTCOLORS(dataSetIndex,:));
        if SANITYCHECKSYMMETRY    
            plot(myLengthSumNewborns{dataSetIndex},1-Ratios{dataSetIndex},'x', 'Color', [.5 .5 .5],'LineWidth',2,'MarkerSize',15);
        end
        
        % cosmetics
        ylim([0,1]);
        xlim([LEFTX,RIGHTX]);

        xlabel(['Parent cell size by ' LENGTHFIELD],'Interpreter','None');
        ylabel('L_{child}/L_{parent}');

        MW_makeplotlookbetter(15);

    end

    disp('section done');

end

%% Now calculate overall histogram
if any(strcmp(RUNSECTIONSFILADIV,{'all','histogramAllLengths'}))

    h=figure(2); clf; hold on;

    % calculate hist
    [count,bincenters] = hist([Ratios{:}],HISTNRBINS);

    % helping lines
    N=5; highestcount=max(count);
    for binIdx=1:N
        for j = 1:(binIdx*2-1)
            plot([(j)/(2*binIdx) (j)/(2*binIdx)],[0, highestcount],'-','Color',[.5 .5 .5],'LineWidth',N-binIdx+1)
        end
    end

    % plot hist (calculated above)
    if exist('histSkel','var') && exist('histArea','var')
        % plot histograms from multiple fields
        l1=plot(histArea(2,:),histArea(1,:),'-','LineWidth',3);    
        l2=plot(histSkel(2,:),histSkel(1,:),'-','LineWidth',3);
        legend([l1,l2],{'Skeleton','Area'});
    else
        % plotting of one histogram
        plot(bincenters,count,'-','LineWidth',3);
    end

    % additional ticks
    plot([0:.1:1],zeros(1,11),'+','Color','k');%,'MarkerFaceColor','k');

    % cosmetics
    %ax=gca; ax.XTick = [0:.1:1];
    ylabel('Count');
    xlabel('L_d/L_p');

    set(gca,'XTick',[0:.1:1])

    ylim([0 highestcount]);

    MW_makeplotlookbetter(15);

    % Save the histogram
    if strcmp(LENGTHFIELD, 'areaPixels')
        histArea=[count;bincenters];
    elseif strcmp(LENGTHFIELD, 'length_skeleton')
        histSkel=[count;bincenters];
    end
    
end

%% Sanity check
if any(strcmp(RUNSECTIONSFILADIV,{'all','sanityParentVsDaughterSum'}))

    % load(datasetsPaths{2})
    figure(3); clf; hold on

    % x=y line
    plot([1,10^6],[1,10^6],'-k')

    % data
    for binIdx=1:numel(datasetsPaths)
        plot(myLengthSumNewborns{dataSetIndex}, myLengthParents{dataSetIndex},'x')
    end

    % cosmetics
    axis equal;
    xlim([0 max([myLengthSumNewborns{:}])]);
    ylim([0 max([myLengthParents{:}])]);
    xlabel('Summed length newborns');
    ylabel('Length of parent');

end
%% Plot historgrams per window
if any(strcmp(RUNSECTIONSFILADIV,{'all','RutgerPlotCombinedHistogram'}))

    NRREGIMES = 4;
    % Set at beginning of script
        %WINDOWBORDERS = [5,10,17,20,30];
        %WINDOWBORDERS = [2,9,16,23,30];
        %WINDOWBORDERS = [2:7:30];
        %WINDOWBORDERS = [3:6:30]; 
    regionsDouble = {[0,0.5],[0,0.5],[0,2/6,.5],[0,2/8,.5],[0,2/10,4/10,.5]}; % boundaries around matching ratios
    LINEWIDTH =2;

    % figure stuff
    h=figure(4); clf; hold on
    myColors = linspecer(numel(WINDOWBORDERS)-1);

    % calculate params
    mybins = linspace(0,1,HISTNRBINS);
    centers = mybins(2:end)-(mybins(2:end)-mybins(1:end-1))/2;
    histData = struct;
    histData.centers = centers; % all datasets and windows share same centers
    for dataSetIndex = 1:numel(datasetsPaths)
        %dataSetIndex=1; % TEMP REMOVE   
        
         % calculated parameters
        if ~exist('RIGHTX','var')
            if strcmp(LENGTHFIELD,'length_skeleton')
                RIGHTX = 30;        
            elseif strcmp(LENGTHFIELD,'areaPixels')
                RIGHTX = 15000;
            else
                RIGHTX  = max(myLengthParents{dataSetIndex})*1.1;
            end
        end

        for windowIndex = 1:(numel(WINDOWBORDERS)-1)

            windowLeft = WINDOWBORDERS(windowIndex);
            windowRight  = WINDOWBORDERS(windowIndex+1);
            windowSize = windowRight-windowLeft;

            % get data
            currentWindowIndices = ...
                find((myLengthSumNewborns{dataSetIndex}>windowLeft) == (myLengthSumNewborns{dataSetIndex}<windowRight));
            currentTotalLength = myLengthSumNewborns{dataSetIndex}(currentWindowIndices);
            currentRatios = Ratios{dataSetIndex}(currentWindowIndices);

            % top row
            % ===
            subplot(2,1,1); hold on
            plot(currentTotalLength,currentRatios,'o','Color',myColors(windowIndex,:),'MarkerFaceColor',myColors(windowIndex,:))
            % cosmetics
            ylim([0,1]);
            xlim([LEFTX,RIGHTX]);
            xlabel(['Parent cell size by ' LENGTHFIELD],'Interpreter','None');
            ylabel('L_{child}/L_{parent}');
            MW_makeplotlookbetter(15)
            title([WHATDATA ' condition']);

            % bottom row with hists
            % ===
            subplot(2,1,2); hold on;    
            % predicted location        
            nrLocations = windowIndex;
            for j = 1:2:(nrLocations*2-1)
                plot([windowLeft, windowRight], [(j)/(2*nrLocations) (j)/(2*nrLocations)],'-','Color',[.5 .5 .5],'LineWidth',2);
            end        
            % histograms
            [counts,edges]=histcounts(currentRatios,mybins); % note histcounts() takes edges as input, hist() does not
            modcounts = (windowSize*counts/sum(counts)+windowLeft);
            plot(modcounts,histData.centers,'-','Color',myColors(windowIndex,:),'LineWidth',LINEWIDTH)
            % cosmetics
            ylim([0,1]);
            xlim([LEFTX,RIGHTX]);
            xlabel(['Histograms (normalized)']);
            ylabel('L_{child}/L_{parent}');
            MW_makeplotlookbetter(15)
            set(gca,'xtick',[])

            % also create stats how many divisions happened at a certain ratio
            eventsPerRatio = histcounts(currentRatios,regionsDouble{windowIndex}); % note histogram() and histcounts() take edges as input, hist() does not

            histData.counts{dataSetIndex}{windowIndex} = counts;
            histData.eventsPerRatio{dataSetIndex}{windowIndex} = eventsPerRatio;

        end
    end

    % Store data and also create summary (mean, sum, normalized) params and store those
    for windowIndex = 1:(numel(WINDOWBORDERS)-1)
        % determine average function
        countsFromMultipleDatasets = arrayfun(@(x) histData.counts{x}{windowIndex}, 1:numel(histData.counts), 'UniformOutput', false);
        histData.meanCounts{windowIndex}=mean(cell2mat(countsFromMultipleDatasets'));
        histData.sumCounts{windowIndex}=sum(cell2mat(countsFromMultipleDatasets'));

        % normalize the pdf
        dt=histData.centers(2)-histData.centers(1);
        histData.normalizedPdf{windowIndex} = histData.sumCounts{windowIndex}./sum(histData.sumCounts{windowIndex})*dt;

        % sum the different event counts for the ratios
        ratiocountsFromMultipleDatasets = arrayfun(@(x) histData.eventsPerRatio{x}{windowIndex}, 1:numel(histData.eventsPerRatio), 'UniformOutput', false);
        histData.sumRatioCounts{windowIndex}=sum(cell2mat(ratiocountsFromMultipleDatasets'));
    end

    REGIONIDX=3;
    dataCounts=cell2mat(arrayfun(@(x) histData.eventsPerRatio{x}{REGIONIDX},1:numel(histData.eventsPerRatio),'UniformOutput',0)');
    myStats = sum(dataCounts)
    disp(['For the ' num2str(REGIONIDX) 'rd region, N is respectively ' num2str(myStats)]);
    
    if ~NOSAVEPLEASE
        saveas(h, [PLOTSAVEDIR WHATDATA '_ratios_histograms.fig']);
        saveas(h, [PLOTSAVEDIR WHATDATA '_ratios_histograms.tif']);
    end        

end
%% Plot simpler figure w. ratios
if any(strcmp(RUNSECTIONSFILADIV,{'all','rutgerPlotv2'}))

    disp(['Plotting N=' num2str(numel([myLengthSumNewborns{:}])) ' independent datapoints.']);
    
    TYPICALDIVSIZE = 7.5; %10; 12;

    % note fig 5 is next one
    h=figure(6); clf; hold on;

    myXlim=[0,max([myLengthSumNewborns{:}])*1.05];
        
    % line of typical size
    xvalues=TYPICALDIVSIZE:.1:myXlim(2);
    plot(xvalues,.5*TYPICALDIVSIZE./xvalues,'--k','LineWidth',2)

    % lines at (1/n)th
    for nn=2:7
        plot(myXlim,[1/nn,1/nn],'-','Color',[.5 .5 .5],'LineWidth',2);
    end

    %{
    % lines at (1/2^n)th
    for nn=2:5
        plot(myXlim,[1/2.^nn,1/2.^nn],'-','Color',[.5 .5 .5],'LineWidth',2);
    end
    %}

    % data
    for dataSetIndex = 1:numel(datasetsPaths)
        plot(myLengthSumNewborns{dataSetIndex},Ratios{dataSetIndex},'o', 'Color', PLOTCOLORS(dataSetIndex,:),'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',PLOTCOLORS(dataSetIndex,:));
        if USESYMMETRY
            plot(myLengthSumNewborns{dataSetIndex},1-Ratios{dataSetIndex},'o', 'Color', PLOTCOLORS(dataSetIndex,:),'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',PLOTCOLORS(dataSetIndex,:));    
        end
    end

    %xlim([LEFTX,RIGHTX]);

    xlabel(['Parent cell size by ' LENGTHFIELD],'Interpreter','None');
    ylabel('L_{child}/L_{parent}');

    xlim(myXlim);
    ylim([0,1]);
    MW_makeplotlookbetter(20);

    if ~NOSAVEPLEASE
        saveas(h, [PLOTSAVEDIR WHATDATA 'ratios_simple.fig']);
        saveas(h, [PLOTSAVEDIR WHATDATA 'ratios_simple.tif']);
        saveas(h, [PLOTSAVEDIR WHATDATA 'ratios_simple.eps'],'epsc');
    end
    
end

%% Plot lifetime against birth length
if any(strcmp(RUNSECTIONSFILADIV,{'all','lifeTimeVsBirthLength1'}))

    ALLDATALOOKSAME=1;

    h=figure(5); clf; hold on
    myColors = linspecer(numel(datasetsPaths));
    myPlotMarkers = 'os^vd<>os^vd<>os^vd<>';

    title([WHATDATA ' condition']);

    %
    meanInterDivisionTimes=[];
    for dataSetIndex = 1:numel(datasetsPaths)

        l=plot(allLengths{dataSetIndex},myLifeTimesSchnitzes{dataSetIndex},'.');    
        set(l,...  
            'Marker',myPlotMarkers(dataSetIndex),...
            'LineWidth',2,...
            'Color',myColors(dataSetIndex,:));%,'MarkerFaceColor',myColors(dataSetIndex,:));


        currentInterDivTimes = myLifeTimesSchnitzes{dataSetIndex};
        meanInterDivisionTimes(dataSetIndex)=mean(currentInterDivTimes(~isnan(currentInterDivTimes)));            
            % note meanInterDivisionTimes is subject to how long cells where in
            % filamented state in this dataset.

    end
    %title(num2str(meanInterDivisionTimes))

    % Plot the means per window (which was set in earlier plot)
    for dataSetIndex = 1:numel(datasetsPaths)

        %allSetsLifeTimes = [lifeTimes{:}];
        %allSetsallLengths = [allLengths{:}];

        meansLengths=[]; meansLifeTimes=[]; stdLifeTimes = [];
        for windowIndex = 1:(numel(WINDOWBORDERS)-1)

            windowLeft = WINDOWBORDERS(windowIndex);
            windowRight  = WINDOWBORDERS(windowIndex+1);
            windowSize = windowRight-windowLeft;

            % get data
            currentWindowIndices = ...
                find((allLengths{dataSetIndex}>windowLeft) == (allLengths{dataSetIndex}<windowRight));
            currentLengths = allLengths{dataSetIndex}(currentWindowIndices);
            currentLifeTimes = myLifeTimesSchnitzes{dataSetIndex}(currentWindowIndices);

            % determine nan values
            indicesBothNoNaN = ~isnan(currentLengths) & ~isnan(currentLifeTimes);

            % determine means
            meansLengths(end+1) = mean(currentLengths(indicesBothNoNaN));
            meansLifeTimes(end+1) = mean(currentLifeTimes(indicesBothNoNaN));   
            stdLifeTimes(end+1) =  std(currentLifeTimes(indicesBothNoNaN));
            stdLifeTimesUncorr(end+1) = stupidstd(currentLifeTimes(indicesBothNoNaN));

        end

        % plot
        errorbar(meansLengths,meansLifeTimes,stdLifeTimes,'ok','MarkerFaceColor','k')
        plot(meansLengths,meansLifeTimes,'-k',...
            'Marker',myPlotMarkers(dataSetIndex),...
            'MarkerFaceColor','k','MarkerSize',10)
    end

    %xlim([LEFTX,RIGHTX]);
    %xlabel(['Birth size by ' LENGTHFIELD ' '],'Interpreter','None');
    xlabel('Birth size (\mum)');
    ylabel('Lifetime [min]');
    MW_makeplotlookbetter(15)

    if ~NOSAVEPLEASE
        saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_histograms.fig']);
        saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_histograms.tif']);
    end

end
%% Similar plot for paper without datasets separate color etc
if any(strcmp(RUNSECTIONSFILADIV,{'all','lifeTimeVsBirthLength2'}))
    hLenLifetime = figure(14); clf; hold on;

    SUBDIR = 'birthsizeLifetime\';

    xdata=[allLengths{:}];
    ydata=[myLifeTimesSchnitzes{:}];
    correctedbdata=[];
    for binIdx=1:numel(datasetsPaths)         
        correctedbdata=[correctedbdata birthTimes{binIdx}-switchTimes(binIdx)];
    end

    selectedPointIndices = correctedbdata>0 & ~isnan(xdata)&~isnan(ydata);
    selectedXdata = xdata(selectedPointIndices);
    selectedYdata = ydata(selectedPointIndices);
    
    myBins = linspace(min(xdata),max(xdata),20);

    %selectedIdxs=~isnan(xdata)&~isnan(ydata);
    %selectedDataX={xdata(selectedIdxs)};
    %selectedDataY={ydata(selectedIdxs)};
    [meanValuesForBinsDynData, binCentersDynData,stdValuesForBinsDynData,stdErrValuesForBins]=binnedaveraging({selectedXdata},{selectedYdata},myBins);

    markerWidth=7; markerColor = [.7 .7 .7];
    hS=scatter(selectedXdata,selectedYdata,...
            7^2,[.5 .5 .5],'filled','MarkerFaceAlpha',.3);
            %'filled','MarkerFaceColor',markerColor,'MarkerEdgeColor','none','MarkerFaceAlpha',1);            
    set(hS, 'SizeData', markerWidth^2);

    if exist('CUSTOMCOLOR')
        theColor=CUSTOMCOLOR;
    else
        theColor='k';
    end
    
    pointsWithMultipleDataIndices = stdValuesForBinsDynData>0;
    errorbar(binCentersDynData(pointsWithMultipleDataIndices),meanValuesForBinsDynData(pointsWithMultipleDataIndices),stdErrValuesForBins(pointsWithMultipleDataIndices),...
        '-o','Color',theColor,'MarkerSize',5,'LineWidth',2);
        %'ok','LineWidth',3,'MarkerFaceColor','k');

    MW_makeplotlookbetter(20);
    xlabel('Birth size (µm)');
    ylabel('Interdivision time (min)');

    if exist('YLIMbirthlife','var')
        ylim(YLIMbirthlife);
    else
        ylim([0,max(ydata)*1.1])
        disp('Set YLIMbirthlife for custom lims');
    end

    if exist('MYLLIM','var')
        xlim(MYLLIM);
    end

    if ~NOSAVEPLEASE
        if ~exist([PLOTSAVEDIR SUBDIR],'dir')
            mkdir([PLOTSAVEDIR SUBDIR]);
        end

        if exist('ONESTOTAKE')
            suffix=num2str(ONESTOTAKE);
            suffix=strrep(suffix,' ','_');
        else
            suffix='';
        end

        saveas(hLenLifetime, [PLOTSAVEDIR SUBDIR WHATDATA '_sizeLifetime_dynamics' suffix '.fig']);
        saveas(hLenLifetime, [PLOTSAVEDIR SUBDIR WHATDATA '_sizeLifetime_dynamics' suffix '.svg']);
        saveas(hLenLifetime, [PLOTSAVEDIR SUBDIR WHATDATA '_sizeLifetime_dynamics' suffix '.png']);
    end
    
    hFig2BD = hLenLifetime;
    
    %%
    figure(104); clf; hold on;
    largeCellsselectedYdata = selectedYdata(selectedXdata>15);
    edges=linspace(0,100,11);
    [n,edges]=histcounts(largeCellsselectedYdata,edges);
    centers=edges(2:end)-(edges(2)-edges(1))/2;
    bar(centers, n);
    xlabel('Interdivision time'); ylabel('Count');
    eventsBigger10=sum(largeCellsselectedYdata>10);
    eventsSmaller10=sum(largeCellsselectedYdata<=10);
    percentageBigger10=100*eventsBigger10./(eventsBigger10+eventsSmaller10);
    disp(['For ' WHATDATA ':']);
    disp(['Of ' num2str(numel(largeCellsselectedYdata)) ' events, ' num2str(percentageBigger10) '% or ' num2str(eventsBigger10) ' events were >10 mins in cells >10µm.']);
    MW_makeplotlookbetter(20);
    
end

%% This plots distributions for each of the bins of the previous plots
% Run previous section first
if any(strcmp(RUNSECTIONSFILADIV,{'all','lifeTimeVsBirthLengthPDFs'}))    
    
    %% Now also get distributions
    figure(101); clf; hold on;
    hFig2Bsub=figure(102); clf; hold on;
    nrYbins=10;
    yBins=linspace(0,100,nrYbins);
    dYbin = (yBins(2)-yBins(1));
    yBinsCenters = yBins(1:end-1)+dYbin/2;
    [pdfValuesForBinsDynData, binCentersDynData]=binnedpdfs({selectedXdata},{selectedYdata},myBins,yBins);
    %pdfsForBins=NaN(numel(edgesForX)-1,numel(edgesForY)-1);
    baseColors = linspecer(3);
    myColors = makeColorMap(baseColors(1,:), baseColors(2,:), baseColors(3,:), (numel(myBins)-1));
    %myColors = makeColorMap([0 1 0], [0 0 1], [0 0 0], (numel(myBins)-1));
    lines = []; legendTitles = {};
    normalizedyBins = {};
    dYNormalizedBin = {};
    normalizedPdf = {};
    normalizedPdfNormBins = {};
    for binIdx = 1:(numel(myBins)-1)
        
        normalizedyBins{binIdx} = yBinsCenters./meanValuesForBinsDynData(binIdx);
        dYNormalizedBin{binIdx} = normalizedyBins{binIdx}(2)-normalizedyBins{binIdx}(1);        
        normalizedPdf{binIdx} = pdfValuesForBinsDynData(binIdx,:)./(sum(pdfValuesForBinsDynData(binIdx,:))*dYbin);
        normalizedPdfNormBins{binIdx} = pdfValuesForBinsDynData(binIdx,:)./(sum(pdfValuesForBinsDynData(binIdx,:))*dYNormalizedBin{binIdx});
        totalCount = sum(pdfValuesForBinsDynData(binIdx,:));
        
        if totalCount>4
            figure(101);
            l=plot(yBinsCenters,normalizedPdf{binIdx},'-','LineWidth',2);            
            %[edges, doubley] = centerstoedges(yBinsCenters,normalizedPdf);
            %l=plot(edges,doubley,'-','LineWidth',2);                        
            set(l,'Color',myColors(binIdx,:));
            
            figure(102);
            l=plot(normalizedyBins{binIdx},normalizedPdfNormBins{binIdx},'-o','LineWidth',2);       
            %[edges, doubley] = centerstoedges(normalizedyBins,normalizedPdf);
            %l=plot(edges,doubley,'-','LineWidth',2);                        
            set(l,'Color',myColors(binIdx,:));

            lines(end+1) = l;
            legendTitles{end+1} = ['' sprintf('%0.1f', binCentersDynData(binIdx)) '?m (N=' num2str(totalCount) ')'];
        end
        
    end
    
    % Cosmetics fig. 101 
    figure(101);
    legend(gca,lines,legendTitles)
    MW_makeplotlookbetter(20);
    xlabel('Interdivision time');    
    ylabel('Normalized probability density');
    % Cosmetics fig. 102 (normalized time, collapsed distr.)
    figure(102);
    legend(gca,lines,legendTitles)
    MW_makeplotlookbetter(20);
    xlabel('Normalized interdivision time');
    ylabel('Normalized probability density');
    %xlim([0,5]);
    
end

%% Plot added length vs. cell length vs. lifetime
if any(strcmp(RUNSECTIONSFILADIV,{'all','adderLengthLifetime'}))
    
    MYXEDGES=linspace(0,20,10)
    MYXEDGES=linspace(0,20,50) % bins for adder plot
    
    %%
    hFig2D=figure(103); clf; hold on;
   
    markerColor=[.5 .5 .5];

    alldatax=[];alldatay=[];
    for datasetIdx=1:numel(datasetsPaths)

        bornAfterSwitchIdxs = (birthTimes{datasetIdx}-switchTimes(datasetIdx))>0;    
        selectedLifeTimeData = myLifeTimesSchnitzes{datasetIdx}(bornAfterSwitchIdxs);

        baseColors = linspecer(3);
        myColors = makeColorMap(baseColors(1,:), baseColors(2,:), baseColors(3,:));
        %{
        lifeTimesNormalized=round(...
            (selectedLifeTimeData-min(selectedLifeTimeData))./    ...
            (max(selectedLifeTimeData)-min(selectedLifeTimeData)) ...
            .*99 ...
            )+1;
        %}
        % Normalize without subtracting min. time
        lifeTimesNormalized=round(...
            (selectedLifeTimeData)./    ...
            (max(selectedLifeTimeData)) ...
            .*99 ...
            )+1;
        lifeTimesNormalized(isnan(lifeTimesNormalized))=1;
        markerColors=myColors(lifeTimesNormalized,:);
        
        l=scatter(  myCellLengthsPerSchnitz{datasetIdx}(bornAfterSwitchIdxs),...
                    myAddedLengthPerSchnitz{datasetIdx}(bornAfterSwitchIdxs),...
                    [],markerColors,'filled','MarkerFaceAlpha',.5);

        %l=scatter(  myCellLengthsPerSchnitz{datasetIdx}(bornAfterSwitchIdxs),...
        %    myAddedLengthPerSchnitz{datasetIdx}(bornAfterSwitchIdxs));
        %set(l,'MarkerFaceColor',markerColors,'MarkerEdgeColor','none','MarkerFaceAlpha',.5);

        alldatax = [alldatax myCellLengthsPerSchnitz{datasetIdx}(bornAfterSwitchIdxs)];
        alldatay = [alldatay myAddedLengthPerSchnitz{datasetIdx}(bornAfterSwitchIdxs)];
        
    end
    
    [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=binnedaveraging({alldatax},{alldatay},MYXEDGES);
    errorbar(binCenters,meanValuesForBins,stdErrValuesForBins,'ok','LineWidth',2);
    
    %
    ylim([0,4]); 
    xlim([0,20]);
    xlabel('Cell length (µm)'); ylabel('Added length \Delta (µm)');
    
    colormap(myColors);
    hColorbar=colorbar;
    
    inputSettings.rangeOut = [0,max(selectedLifeTimeData)];
    inputSettings.desiredSpacing = 25;
    inputSettings.rangeIn = [0,1];
    inputSettings.desiredDecimalsTicks=0;
    [tickLocationsOldNewMetric, correspdongingLabels] = labelremapping(inputSettings);
    set(hColorbar,'YTick',tickLocationsOldNewMetric,'YTickLabel',correspdongingLabels);
    ylabel(hColorbar,'Interdivision time (mins)');

    MW_makeplotlookbetter(20);
    
    %% Now also create pdfs
    MYXEDGES=MYXEDGES;
    xCenters = MYXEDGES(1:end-1)+(MYXEDGES(2)-MYXEDGES(1))/2;
    yEdges=linspace(0,4,10);
    dY = (yEdges(2)-yEdges(1));
    yCenters = yEdges(1:end-1)+dY/2;
    [pdfsForBins, binCenters]=binnedpdfs({alldatax},{alldatay},MYXEDGES,yEdges);
    
    figure; clf; hold on;
    legendTitles={}; lHandles=[];
    for binIdx=1:size(pdfsForBins,1)
        currentPdf = pdfsForBins(binIdx,:);
        totalN=sum(currentPdf)
        
        dYNorm          = dY./meanValuesForBins(binIdx);
        yCentersNorm    = yCenters./meanValuesForBins(binIdx);
        
        currentPdfNorm = currentPdf./(totalN.*dYNorm);
        
        if totalN>10
            lHandles(end+1)=plot(yCentersNorm,currentPdfNorm,'LineWidth',2);
        end
        
        legendTitles{end+1} = ['L_b=' sprintf('%0.1f', xCenters(binIdx)) '?m (N=' num2str(totalN) ')'];
        
    end
    
    legend(lHandles,legendTitles);
    xlabel('\Delta/<\Delta>');
    ylabel('Probability');
    MW_makeplotlookbetter(20);
    
end

%% Plot with time-coded ratios
if any(strcmp(RUNSECTIONSFILADIV,{'all','rutgerPlotFinal'}))

    % Note: there was once also a 2nd figure to plot a selection, but this is
    % now covered by a later figure (with color coding) and has been commented
    % out here..

    %%
    ALPHA = .15; %.15;
    
    MINUTELIFETIMETRESHOLD=20;
    if ~exist('PLOTHELPINGLINES','var')
        PLOTHELPINGLINES=1;
    end
    if ~exist('COLORCODEMARKERS','var')
        COLORCODEMARKERS=1; 
    end

    h1=figure(7); clf; hold on
    %h2=figure(8); clf; hold on

    % Get a colormap (is 64 long by default)
    timeColormap = colormap(jet); % winter
    rowsInColormap = size(timeColormap,1);

    % 
    longestTime = max([myLifeTimesSchnitzes{:}]);
    %longestTime = max([birthTimes{:}]);
    fromLifeTimeToColorCodeFactor = longestTime/(rowsInColormap-1);

    divisionNrCounter=0;
    for dataSetIndex = 1:numel(datasetsPaths)

        for binIdx = 1:numel(myLengthSumNewborns{dataSetIndex})

            %currentLifeTime = lifeTimes{dataSetIndex}(myNewBornSchnitzNrs{dataSetIndex}(i)); % daughter lifetime
            currentTime = myLifeTimeParents{dataSetIndex}(binIdx); % parent lifetime
            %currentTime = birthTimes{dataSetIndex}(i); % time at which born

            % determine color if lifetime known
            if isnan(currentTime) || ~COLORCODEMARKERS
                colorForThisDataPoint = [.7 .7 .7];
                %colorForThisDataPoint = [.2 .2 .2];
            else    
                colorForThisDataPoint = ...
                    timeColormap(...
                        round(currentTime/fromLifeTimeToColorCodeFactor)+1 ...
                        ,:);
            end

            % plot all datapoint to one figure
            figure(h1.Number);
            scatter(myLengthSumNewborns{dataSetIndex}(binIdx),Ratios{dataSetIndex}(binIdx),12^2,'filled',...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',colorForThisDataPoint,...
                    'MarkerFaceAlpha',ALPHA,...
                    'MarkerEdgeAlpha',ALPHA);
            if USESYMMETRY
                scatter(myLengthSumNewborns{dataSetIndex}(binIdx),1-Ratios{dataSetIndex}(binIdx),12^2,'filled',...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',colorForThisDataPoint,...
                    'MarkerFaceAlpha',ALPHA,...
                    'MarkerEdgeAlpha',ALPHA);
            end
    
            divisionNrCounter=divisionNrCounter+numel(myLengthSumNewborns{dataSetIndex}(binIdx));

            % plot selection of datapoints to 2nd figure
            %{
            figure(h2.Number);
            if currentTime>MINUTELIFETIMETRESHOLD % || isnan(currentLifeTime)
                plot(myLengthSumNewborns{dataSetIndex}(i),Ratios{dataSetIndex}(i),'o','MarkerSize',10,...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',colorForThisDataPoint)
                if USESYMMETRY
                    plot(myLengthSumNewborns{dataSetIndex}(i),1-Ratios{dataSetIndex}(i),'o','MarkerSize',10,...
                        'MarkerEdgeColor',colorForThisDataPoint,...
                        'MarkerFaceColor',colorForThisDataPoint)            
                end
            end
            %}
        end

    end

    for myfignum = h1.Number%:h2.Number
        figure(myfignum);

        xlabel('Length of mother cell Lm (µm)');
        ylabel('Relative division location S');
        
        %xlabel('Length of mother cell \it{L_m}\rm (?m)');
        %ylabel('Relative division location \itS\rm');

        MW_makeplotlookbetter(20);

        if COLORCODEMARKERS
            hc = colorbar;
            %ax = gca;
            %ax.XTickLabel = {'1000','1'};
            caxis([0,longestTime])
            ylabel(hc, 'Daughter interdivision time')
        end

    end

    if exist('MYXLIM','var')
        figure(7); xlim(MYXLIM); %figure(8); xlim(MYXLIM);
    else
        disp('Set MYXLIM to [limleft limright] to adjust plots x axis to custom values');
    end

    % plot helping lines at 1/2n
    if PLOTHELPINGLINES
        %%
        N=5;
        for windowIndex=1:(numel(WINDOWBORDERS)-1)
            for j = 1:2:(windowIndex*2-1)
                if ~exist('THINHELPINGLINES')
                    plot([WINDOWBORDERS(windowIndex), WINDOWBORDERS(windowIndex+1)], [(j)/(2*windowIndex) (j)/(2*windowIndex)],'-','Color',BARCOLORS(windowIndex,:),'LineWidth',4)
                else
                    plot([WINDOWBORDERS(windowIndex), WINDOWBORDERS(windowIndex+1)], [(j)/(2*windowIndex) (j)/(2*windowIndex)],':','Color',BARCOLORS(windowIndex,:),'LineWidth',4)
                end
            end
        end
    end        
    
    if ~exist('NOBIRTHLINE') 
        if exist('MYXLIM')
            %DAUGHTERSIZE=4;
            for DAUGHTERSIZE=[2,4]
                xvalues=linspace(2*DAUGHTERSIZE,MYXLIM(2),100);
                l=plot(xvalues,DAUGHTERSIZE./xvalues);
                set(l, 'Color','k','LineStyle','--','LineWidth',3);
            end
            
            ylim([0,1]);
        else
                warning('Not plotting birthline since MYXLIM not set.');
        end
        
    end

    if ~NOSAVEPLEASE
        MW_makeplotlookbetter(12*2);

        SIZE=[6.8 6.8];

        h1.Units = 'centimeters';
        h1.Position = [0 0 SIZE]*2+1;    

        h1.PaperUnits = 'centimeters';
        h1.PaperPosition = [0 0 SIZE]*2;

        saveas(h1, [PLOTSAVEDIR WHATDATA '_ratiosByInterdiv.fig']);
        saveas(h1, [PLOTSAVEDIR WHATDATA '_ratiosByInterdiv.tif']);

        % For some reason Matlab refuses to produce a vector image here?!
        saveas(h1, [PLOTSAVEDIR WHATDATA '_ratiosByInterdiv.svg']);
        saveas(h1, [PLOTSAVEDIR WHATDATA '_ratiosByInterdiv.eps'],'epsc');
        saveas(h1, [PLOTSAVEDIR WHATDATA '_ratiosByInterdiv.pdf']);

        %saveas(h2, [PLOTSAVEDIR WHATDATA '_ratiosByInterdivSelected.fig']);
        %saveas(h2, [PLOTSAVEDIR WHATDATA '_ratiosByInterdivSelected.tif']);
    end
    
    figure1DHandle = h1;
    hSupp5 = h1;
    
    disp(['N=' num2str(divisionNrCounter) ' divisions plotted.']);
    
end
%% Histogram of lifetimes
if any(strcmp(RUNSECTIONSFILADIV,{'all','histogramLifetimes'}))

    h=figure(9); clf; hold on

    hist([myLifeTimeParents{:}],30)
    xlabel('Interdivision lifetime [min]');
    ylabel('Count')

    MW_makeplotlookbetter(20);

    if ~NOSAVEPLEASE
        saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_generalhistogram.fig']);
        saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_generalhistogram.tif']);
    end
    
end

%% histogram of cell sizes
if any(strcmp(RUNSECTIONSFILADIV,{'all','histogramCellSizes'}))
    
    figure(10); clf; hold on;
    allLengthsPile = [schnitzcells(:).(LENGTHFIELD)];
    histogram(allLengthsPile,100);

    xlabel('Length (um)');
    ylabel('Count');
    title('Cell lenghts observed during movie');
    MW_makeplotlookbetter(15);
    
end

%% Plot with ratios, black/white time coding
if any(strcmp(RUNSECTIONSFILADIV,{'all','rutgerPlotBlackWhiteTimecoding'}))

    disp(['Plotting N=' num2str(numel([myLengthSumNewborns{:}])) ' independent datapoints.']);
    
    MINUTELIFETIMETRESHOLD=20;

    h1=figure(11); clf; hold on

    % Get a colormap (is 64 long by default)
    timeColormap = colormap(jet); % winter
    rowsInColormap = size(timeColormap,1);

    % 
    longestTime = max([myLifeTimesSchnitzes{:}]);
    %longestTime = max([birthTimes{:}]);
    fromLifeTimeToColorCodeFactor = longestTime/(rowsInColormap-1);

    for dataSetIndex = 1:numel(datasetsPaths)

        for cellIdx = 1:numel(myLengthSumNewborns{dataSetIndex})            

            %currentLifeTime = lifeTimes{dataSetIndex}(myNewBornSchnitzNrs{dataSetIndex}(i)); % daughter lifetime
            currentTime = myLifeTimeParents{dataSetIndex}(cellIdx); % parent lifetime
            %currentTime = birthTimes{dataSetIndex}(i); % time at which born

            % determine color if lifetime known
            if isnan(currentTime)
                colorForThisDataPoint = [0 0 0];
                faceForPoint = 'None';
            else    
                if currentTime<MINUTELIFETIMETRESHOLD
                    colorForThisDataPoint = [0 0 0];
                    faceForPoint = 'None';
                else
                    colorForThisDataPoint = [0 0 0];            
                    faceForPoint = 'k';
                end
            end

            % plot all datapoint to one figure
            plot(myLengthSumNewborns{dataSetIndex}(cellIdx),Ratios{dataSetIndex}(cellIdx),'o','MarkerSize',2,...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',faceForPoint)
            if USESYMMETRY
                plot(myLengthSumNewborns{dataSetIndex}(cellIdx),1-Ratios{dataSetIndex}(cellIdx),'o','MarkerSize',2,...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',faceForPoint)
            end

        end

    end

    % plot helping lines at 1/2n
    if PLOTHELPINGLINES
        N=5;
        for windowIndex=1:(numel(WINDOWBORDERS)-1)
            for j = 1:2:(windowIndex*2-1)
                plot([WINDOWBORDERS(windowIndex), WINDOWBORDERS(windowIndex+1)], [(j)/(2*windowIndex) (j)/(2*windowIndex)],'-','Color',BARCOLORS(windowIndex,:),'LineWidth',8)
            end
        end
    end 
    
    %xlabel('Summed daughter length (µm)');
    xlabel('Mother length (µm)');
    %ylabel('Relative division location (S)');
    ylabel('Division location, S');

    MW_makeplotlookbetter(20);

    if exist('MYLLIM','var')
        figure(h1.Number); xlim(MYLLIM);
    else
        disp('Set MYLLIM to [limleft limright] to adjust plots x axis to custom values');
    end

    if ~NOSAVEPLEASE

        if ~exist([PLOTSAVEDIR 'divPlot\'])
            mkdir([PLOTSAVEDIR 'divPlot\']);
        end

        if exist('ONESTOTAKE')
            suffix=num2str(ONESTOTAKE);
            suffix=strrep(suffix,' ','_');
        else
            suffix='';
        end

        saveas(h1, [PLOTSAVEDIR 'divPlot\' WHATDATA  '_newRutgerPlot' suffix '.fig']);
        saveas(h1, [PLOTSAVEDIR 'divPlot\' WHATDATA  '_newRutgerPlot' suffix '.tif']);
        saveas(h1, [PLOTSAVEDIR 'divPlot\' WHATDATA  '_newRutgerPlot' suffix '.png']);
    end

    hSupp5 = h1;
    
end
%% Size time evolution plot
if any(strcmp(RUNSECTIONSFILADIV,{'all','sizeTraces'}))

    %%
    hSizeTime = figure(12); clf; hold on;

    MARKERSIZE=8;
    LINEWIDTH=6;
    FONTSIZE=20;
    SUBDIR = 'trace_dynamics\';
    STARTSCHNITZES = 1;
    COLORDOTS = 0; % 1 or 0
    MAXTRACESHIGHLIGHTED=4;

    traceColors = linspecer(numel(STARTSCHNITZES)*numel(datasetsPaths));

    myTraces={}; myDivisions = {};
    allTimeScatterData={}; allLengthScatterData={};
    endOfLineages ={};
    % loop over datasets
    for dataSetIndex = 1:numel(datasetsPaths)

        % Get the data
        if ~exist('simulatedschnitzcells','var')
            if ~SPECIALCASE
                schnitzcells = loadandrename(datasetsPaths{dataSetIndex});            
            else
                schnitzcells=S_all_shifted{dataSetIndex};
            end
        end

        timeScatterData   = [schnitzcells(:).(TIMEFIELD)];
        lengthScatterData = [schnitzcells(:).(LENGTHFIELD)];
        allTimeScatterData{dataSetIndex} = timeScatterData;
        allLengthScatterData{dataSetIndex} = lengthScatterData;

        % plotting all data from this dataset
        if COLORDOTS
            markerColor = traceColors(dataSetIndex,:); markerColor = 1-(1-markerColor)/3;
        else
            markerColor=[.7 .7 .7]; 
        end
        markerWidth = 4; % 8
        
        
        % plot the lineages
        % Get all traces of the lineages
        [myTracesAll] = MW_gettracefromschnitzcellsreverse(schnitzcells,[],'time',{LENGTHFIELD});
        for idxTrace = 1:numel(myTracesAll)
            l=plot(myTracesAll(idxTrace).time-switchTimes(dataSetIndex),myTracesAll(idxTrace).(LENGTHFIELD),'-','Color',markerColor,'LineWidth',2);
            % change markerColor to [markerColor .1] to set alpha, or use:
            %l.Color(4)=.1;
            % however, multiple traces are plotted on top of each other so
            % it remains hard to see.
        end
        
        
        % plot the dots
        hS=scatter(timeScatterData-switchTimes(dataSetIndex),lengthScatterData,'filled',...
                'MarkerFaceColor',markerColor,'MarkerEdgeColor','none','MarkerFaceAlpha',1);
                %'MarkerFaceColor',markerColor,'MarkerEdgeColor','none','MarkerFaceAlpha',.1); % alpha
        set(hS, 'SizeData', markerWidth^2);

        % Now get ends of schnitzcells lineages
        endSchnitzesIdxs = find([schnitzcells(:).D]==0);
        theEndTimes = arrayfun(@(x) schnitzcells(x).(TIMEFIELD)(end), endSchnitzesIdxs);
        theEndLengths = arrayfun(@(x) schnitzcells(x).(LENGTHFIELD)(end), endSchnitzesIdxs);
        endOfLineages{dataSetIndex}= [theEndTimes; theEndLengths];       
        
        % Now also make example traces
        if strcmp(WHATDATA,'temperature')
            % [myCurrentTraces, myCurrentDivisions,randomSelectedTracesIdxs] = MW_gettracefromschnitzcells(schnitzcells,STARTSCHNITZES,TIMEFIELD,LENGTHFIELD,'random'); 
            % randomSelectedTracesIdxs
            selectedTraces = {[2     5    11    19    82    89   166   241],[6     7    76]};
            [myCurrentTraces, myCurrentDivisions,traceIdxs] = MW_gettracefromschnitzcells(schnitzcells,STARTSCHNITZES,TIMEFIELD,LENGTHFIELD,selectedTraces{dataSetIndex});
        else
            [myCurrentTraces, myCurrentDivisions] = MW_gettracefromschnitzcells(schnitzcells,STARTSCHNITZES,TIMEFIELD,LENGTHFIELD,'longest');
        end
        myTraces{dataSetIndex} = myCurrentTraces;
        myDivisions{dataSetIndex} = myCurrentDivisions;                    

    end

    % Now indicate where a dataset ends
    for dataSetIndex = 1:numel(datasetsPaths)
        plot(endOfLineages{dataSetIndex}(1,:)-switchTimes(dataSetIndex),endOfLineages{dataSetIndex}(2,:),...
                'sk','MarkerSize',MARKERSIZE);
    end

    % for last dataset make additional trace following shortest schnitz
    dataSetIndex=min(1,numel(datasetsPaths));

    % Load schnitzcells param
    if ~SPECIALCASE
        schnitzcells = loadandrename(datasetsPaths{dataSetIndex});            
    else
        schnitzcells=S_all_shifted{dataSetIndex};
    end

    [myShortTrace, myShortTraceDivisions] = MW_gettracefromschnitzcells(schnitzcells,STARTSCHNITZES,TIMEFIELD,LENGTHFIELD,'shortest');
    % now plot traces that were made
    % following short schnitzes
    startSchnitzIdx=1;
    plot(myShortTrace{startSchnitzIdx}(:,1)-switchTimes(dataSetIndex),myShortTrace{startSchnitzIdx}(:,2),'Color','k','LineWidth',LINEWIDTH);
    plot(myShortTraceDivisions{startSchnitzIdx}(:,1)-switchTimes(dataSetIndex),myShortTraceDivisions{startSchnitzIdx}(:,2),'o','Color','k','MarkerSize',MARKERSIZE,'MarkerFaceColor','k');
    % following long schnitzes
    counter=0;
    for dataSetIndex = 1:min(numel(datasetsPaths),MAXTRACESHIGHLIGHTED)
    for startSchnitzIdx = 1:STARTSCHNITZES
        counter=counter+1;
        plot(myTraces{dataSetIndex}{startSchnitzIdx}(:,1)-switchTimes(dataSetIndex),myTraces{dataSetIndex}{startSchnitzIdx}(:,2),'Color',traceColors(counter,:),'LineWidth',LINEWIDTH);
        plot(myDivisions{dataSetIndex}{startSchnitzIdx}(:,1)-switchTimes(dataSetIndex),myDivisions{dataSetIndex}{startSchnitzIdx}(:,2),'o','Color',traceColors(counter,:),'MarkerSize',MARKERSIZE,'MarkerFaceColor',traceColors(counter,:));
    end
    end

    % cosmetics
    MW_makeplotlookbetter(FONTSIZE);
    xlabel('Time (min)');
    ylabel('Cell length (µm)');
    %xlim([0,120]);

    if exist('MYTLIM','var')
        disp('set MYTLIM for fig 12');
        xlim(MYTLIM);
    end

    if ~NOSAVEPLEASE
        if ~exist([PLOTSAVEDIR SUBDIR],'dir')
            mkdir([PLOTSAVEDIR SUBDIR]);
        end

        if exist('ONESTOTAKE')
            suffix=num2str(ONESTOTAKE);
            suffix=strrep(suffix,' ','_');
        else
            suffix='';
        end
        if exist('COLORDOTS')
            suffix=[suffix '_colDots'];
        end

        saveas(hSizeTime, [PLOTSAVEDIR SUBDIR WHATDATA '_colLength_dynamics' suffix '.fig']);
        saveas(hSizeTime, [PLOTSAVEDIR SUBDIR WHATDATA '_colLength_dynamics' suffix '.svg']);
        saveas(hSizeTime, [PLOTSAVEDIR SUBDIR WHATDATA '_colLength_dynamics' suffix '.png']);
    end
    
    hFig2AC = hSizeTime;
    
    % legend('1','2','3','4','5')
    
end
%% 2nd plot for inset
% ========
if any(strcmp(RUNSECTIONSFILADIV,{'all','sizeTraces'})) % part 2 for inset

    hSizeTimeInset = figure(13); clf; hold on;

    SUBDIR = 'trace_dynamics\';
    INSETMARKERWIDTH=14;

    markerColor = [.7 .7 .7];
    lSinset=scatter([allTimeScatterData{:}], [allLengthScatterData{:}],...
                'filled',...
                'MarkerFaceColor',markerColor,'MarkerEdgeColor','none','MarkerFaceAlpha',1);
    set(lSinset, 'SizeData', (INSETMARKERWIDTH)^2);

    % now plot traces that were made
    counter=0;
    for dataSetIndex = 1:numel(datasetsPaths)
    for startSchnitzIdx = 1:STARTSCHNITZES    
        counter=counter+1;
        plot(myTraces{dataSetIndex}{startSchnitzIdx}(:,1),myTraces{dataSetIndex}{startSchnitzIdx}(:,2),'Color',traceColors(counter,:),'LineWidth',6);
        plot(myDivisions{dataSetIndex}{startSchnitzIdx}(:,1),myDivisions{dataSetIndex}{startSchnitzIdx}(:,2),'o','Color',traceColors(counter,:),'MarkerSize',INSETMARKERWIDTH,'MarkerFaceColor',traceColors(counter,:));
    end
    end

    xlim([0,40]);
    ylim([0,max([allLengthScatterData{:}])]);

    % cosmetics
    MW_makeplotlookbetter(30);

    if ~NOSAVEPLEASE
        saveas(hSizeTimeInset, [PLOTSAVEDIR SUBDIR WHATDATA '_colLength_dynamicsINSET.fig']);
        saveas(hSizeTimeInset, [PLOTSAVEDIR SUBDIR WHATDATA '_colLength_dynamicsINSET.svg']);
    end

end

%% Normalize birth size vs. lifetime with # potential division sites
if any(strcmp(RUNSECTIONSFILADIV,{'all','birthSizeLifeTimeRingNorm'}))

    hSupp2e=figure(20); clf; hold on;

    %regimeBounds = [0,7.35,13.59,20];
    regimeBounds = [0,7.35,13.97,20];
    ringCounts = [1,2,3];

    myColorsClusters = linspecer(numel(ringCounts));

    %normalizedByRingsYData = NaN(size(selectedYdata));
    normalizedByRingsYData={}; normalizedByRingsXData={};
    for regimeIdx=1:(numel(regimeBounds)-1)



        currentBounds = regimeBounds(regimeIdx:regimeIdx+1);
        idxToTransfom = selectedXdata>currentBounds(1) & selectedXdata<currentBounds(2);    
        normalizedByRingsYData{regimeIdx} = selectedYdata(idxToTransfom).*ringCounts(regimeIdx);
        normalizedByRingsXData{regimeIdx} = selectedXdata(idxToTransfom);

        scatter(normalizedByRingsXData{regimeIdx},normalizedByRingsYData{regimeIdx},'filled',...
                'MarkerFaceColor',myColorsClusters(regimeIdx,:),'MarkerEdgeColor','none','MarkerFaceAlpha',1);

    end

    %


    for regimeIdx=1:(numel(regimeBounds)-1)
        [meanValuesForBinsDynData, binCentersDynData,stdValuesForBinsDynData,stdErrValuesForBins]=binnedaveraging({normalizedByRingsXData{regimeIdx}},{normalizedByRingsYData{regimeIdx}},myBins);
        errorbar(binCentersDynData,meanValuesForBinsDynData,stdValuesForBinsDynData,'ok-','LineWidth',3,'MarkerFaceColor','k');
    end


    % cosmetics

    xlabel('Birth size (um)');
    ylabel(['Interdivision time (mins)' 10 'multiplied by ring count']);

    MW_makeplotlookbetter(20);

    % Alternatively, divide by times unit size

    figure(21); clf; hold on;

    lenNormDataY = selectedYdata.*sqrt(selectedXdata);
    lenNormDataY = selectedYdata.*selectedXdata;
    scatter(selectedXdata,lenNormDataY,'filled',...
            'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','MarkerFaceAlpha',1);
    [meanValuesForBinsDynData, binCentersDynData,stdValuesForBinsDynData,stdErrValuesForBins]=binnedaveraging({selectedXdata},{lenNormDataY},myBins);
        errorbar(binCentersDynData,meanValuesForBinsDynData,stdValuesForBinsDynData,'ok-','LineWidth',3,'MarkerFaceColor','k');    

    ylim([0,150]);
    xlim([0,20]);    

    xlabel('Birth size (um)');
    %ylabel(['Interdivision time (mins)' 10 'multiplied by sqrt size']);
    ylabel(['Interdivision time (mins)' 10 'multiplied by size']);

    MW_makeplotlookbetter(20);

end    

%% now re-normalize temperature interdivision data by actual ring count
if any(strcmp(RUNSECTIONSFILADIV,{'all','interdivVsLengthNormalizedWithRingCount'}))

    threeColors=linspecer(3);
    
    % determine binning for these plots
    ringBinedges = [0:2:40];
    
    % We want all interdiv. statistics which are put in combinedDynamicsData
    % in script20161221_filarecovery_birthSizeLifeTimeDynamics_all.m
    
    % First select temperature data
    dataXtemp=combinedDynamicsData.(WHATDATA).selectedXdata;
    dataYtemp=combinedDynamicsData.(WHATDATA).selectedYdata;                         
    
    % Create a figure that shows the distribution of rings
    hFigS7 = [];    
    hFigS7(1) = figure; clf; hold on;
    
    scatter(lengthsPerCell,peakCountsPerCell,5^2,[.7 .7 .7],'MarkerFaceColor',[.7 .7 .7]);   
    
    % mean
    [meanValuesForBinsDynData, binCentersDynData,stdValuesForBinsDynData,stdErrValuesForBins]=binnedaveraging({lengthsPerCell},{peakCountsPerCell},ringBinedges);
    plot(binCentersDynData,meanValuesForBinsDynData,'ok','MarkerFaceColor','k','MarkerSize',10);    
    
    % fit
    selectedIndices = lengthsPerCell<40;
    ringNrFit = polyfit(lengthsPerCell(selectedIndices),peakCountsPerCell(selectedIndices),1)
    plot([0:50],[0:50]*ringNrFit(1)+ringNrFit(2),'r-');
    
    xlabel('Length of cell'); ylabel('Number of rings');
    MW_makeplotlookbetter(20);
    
    % Alternative way of getting averages of data above:
    % the following data (peakCountsPerCell, lengthsPerCell) is collected  in 
    % script20160422_filamentRecoveryFtslabelLocations analyzeData
    %{
    ringPdf   = histcounts(peakCountsPerCell,ringBinedges);
    lengthPdf = histcounts(lengthsPerCell,ringBinedges);
    ringCountsPerBacterium = ringCounts./lengthPdf;
    
    plot(ringBins,ringCountsPerBacterium,'ok','MarkerFaceColor','k','MarkerSize',10);
    %}
    
    %{
    % create matching normalization array
    %ringCountsNormalized = ringCountsPerBacterium./max(ringCountsPerBacterium); has inf value unfortunately
    normalizationArray = NaN(size(dataXtemp));
    myBinColors = linspecer(numel(ringBins));
    fillColorArray = NaN(size(dataXtemp));
    for idx = 1:numel(dataXtemp)
        % find which bin it belongs to 
        delta=abs(ringBins-dataXtemp(idx));
        correspondingBinIdx=find(min(delta)==delta);
        correspondingRingCount = ringCountsPerBacterium(correspondingBinIdx);

        normalizationArray(idx) = correspondingRingCount;
        % cosmetic to show which bin was used
        fillColorArray(idx) = myBinColors(correspondingBinIdx);
    end

    normalizeddataYtemp = dataYtemp.*normalizationArray;   
    
    hFigS7C = figure; clf; hold on;
    scatter(dataXtemp,normalizeddataYtemp,7^2,fillColorArray);
    xlim([0,40]); 
    ylim([0,max(normalizeddataYtemp)]); 
    xlabel('Length of cell (?m)');
    ylabel('Interdivision time * ring count');
    MW_makeplotlookbetter(20);

    [meanValuesForBinsDynData, binCentersDynData,stdValuesForBinsDynData,stdErrValuesForBins]=binnedaveraging({dataXtemp},{normalizeddataYtemp},ringBinedges);
    errorbar(binCentersDynData,meanValuesForBinsDynData,stdValuesForBinsDynData,'ok-','LineWidth',3,'MarkerFaceColor','k');    
    
    % 
    figure(); clf; hold on;
    [meanValuesForBinsDynData, binCentersDynData,stdValuesForBinsDynData,stdErrValuesForBins]=binnedaveraging({dataXtemp},{dataYtemp},ringBinedges);    
    errorbar(binCentersDynData,meanValuesForBinsDynData.*binCentersDynData,stdValuesForBinsDynData,'ok-','LineWidth',3,'MarkerFaceColor','k');    
    %}        
    
    %% fit the temperature w. 1/L
    hFigS7(2) = figure; clf; hold on;
    l=plot(dataXtemp, dataYtemp,'o','Color',threeColors(2,:),'MarkerFaceColor',threeColors(2,:),'MarkerSize',3);
    legend(l,'Temperature');
    [Ymean, Xmean, stdValuesForBinsDynData,stdErrValuesForBins]=binnedaveraging({dataXtemp},{dataYtemp},ringBinedges);    
    selectedIdx = dataXtemp>0;% & dataXtemp<20;
    fitPowersTemp = polyfit(dataXtemp(selectedIdx),dataYtemp(selectedIdx).*dataXtemp(selectedIdx),1)
    
    plot(Xmean,Ymean,'ok','MarkerFaceColor','k','MarkerSize',10);
    %plot([0:60],([0:60].*fitPowers(1)+fitPowers(2))./[0:60],'k:','LineWidth',5);
    plot([0:60],fitPowersTemp(1)+fitPowersTemp(2)./[0:60],'k:','LineWidth',5);
    
    xlim([0,40]);
    ylim([0,100]);
    
    xlabel('Length of cell (?m)');
    ylabel('Interdivision time');
    MW_makeplotlookbetter(20);
    
    %% fit the tetracycline data    
    dataXtet=combinedDynamicsData.('tetracycline').selectedXdata;
    dataYtet=combinedDynamicsData.('tetracycline').selectedYdata;
    [tetYmean, tetXmean, stdValuesForBinsDynData,stdErrValuesForBins]=binnedaveraging({dataXtet},{dataYtet},ringBinedges);    
    selectedIdx = dataXtet>0;% & dataXtet<20;
    fitPowersTet = polyfit(dataXtet(selectedIdx),dataYtet(selectedIdx).*dataXtet(selectedIdx),1)
    %{
    % Plot transformed data
    figure; clf; hold on; 
    plot(dataXtet,dataYtet.*dataXtet,'.');
    plot([0:60],[0:60].*fitPowersTet(1)+fitPowersTet(2),'-');
    %}
    
    % plot with 
    hFigS7(3) = figure; clf; hold on;
    l=plot(dataXtet, dataYtet,'o','Color',threeColors(1,:),'MarkerFaceColor',threeColors(1,:),'MarkerSize',3)
    legend(l,'Tetracycline');
    plot(tetXmean,tetYmean,'o','Color','k','MarkerFaceColor','k','MarkerSize',10);
    %plot([0:60],([0:60].*fitPowers(1)+fitPowers(2))./[0:60],'k:','LineWidth',5);
    plot([0:60],fitPowersTet(1)+fitPowersTet(2)./[0:60],'k:','LineWidth',5);
    
    xlim([0,40]);
    ylim([0,100]);
    
    xlabel('Length of cell (?m)');
    ylabel('Interdivision time');
    MW_makeplotlookbetter(20);      
    
    %% fit sulA to 1/L
    hFigS7(4) = figure; clf; hold on;
    dataXsul=combinedDynamicsData.('sulA').selectedXdata;
    dataYsul=combinedDynamicsData.('sulA').selectedYdata;
    
    l=plot(dataXsul, dataYsul,'o','Color',threeColors(3,:),'MarkerFaceColor',threeColors(3,:),'MarkerSize',3);
    legend(l,'SulA');
    [Ymean, Xmean, stdValuesForBinsDynData,stdErrValuesForBins]=binnedaveraging({dataXsul},{dataYsul},ringBinedges);    
    selectedIdx = dataXsul>0;% & dataXsul<20;
    fitPowersTemp = polyfit(dataXsul(selectedIdx),dataYsul(selectedIdx).*dataXsul(selectedIdx),1)
    
    plot(Xmean,Ymean,'ok','MarkerFaceColor','k','MarkerSize',10);
    %plot([0:60],([0:60].*fitPowers(1)+fitPowers(2))./[0:60],'k:','LineWidth',5);
    plot([0:60],fitPowersTemp(1)+fitPowersTemp(2)./[0:60],'k:','LineWidth',5);
    
    xlim([0,40]);
    ylim([0,100]);
    
    xlabel('Length of cell (?m)');
    ylabel('Interdivision time');
    MW_makeplotlookbetter(20);
    
end
%%
if any(strcmp(RUNSECTIONSFILADIV,{'all','SlocationAgainstTime'}))
        
    if ~exist('MYSUBPLOT')
        MYSUBPLOT=1;
    end
    
    hSvsT=figure(5); 
    
    % Create subplot for supplemental figure, depending on flag TIMEDIVPLOTMAINFIG
    if exist('TIMEDIVPLOTMAINFIG','var')
        hSvsT=figure(5); clf; hold on;
    else
        hSvsT=figure(5); 
        if MYSUBPLOT==1, clf; hold on; end;    
        subplot(1,3,MYSUBPLOT); hold on; 
    end
    
    NrDataSets = numel(datasetsPaths);
    
    markerColor=[.5 .5 .5];
    if exist('TIMEDIVPLOTMAINFIG','var')
        markerWidth = 5;
    else
        markerWidth = 10;
    end
    
    ALPHA = .3;
    myColors = linspecer(NrDataSets);   

    allTimes=[]; % just for calculating xlim
    for dataSetIndx = 1:NrDataSets

        % In suppl. fig color data per dataset, don't in main figure.
        if ~exist('TIMEDIVPLOTMAINFIG','var')
            markerColor=myColors(dataSetIndx,:);
        end

        tData = birthTimes{dataSetIndx}(find(usedSchnitzes{dataSetIndx}))-switchTimes(dataSetIndx);
        sData = Ratios{dataSetIndx};

        hS=scatter(tData,sData,'filled',...
                'MarkerFaceColor',markerColor,'MarkerEdgeColor','none','MarkerFaceAlpha',1);
        set(hS, 'SizeData', markerWidth^2,...
                'MarkerFaceAlpha',ALPHA,...
                'MarkerEdgeAlpha',ALPHA);

        allTimes = [allTimes tData];
    end

    if exist('TIMEDIVPLOTMAINFIG','var')
        xlabel('Time (mins)');
        ylabel('Location S');
    else
        if MYSUBPLOT==2
            xlabel('Time [minutes]');
        end
        if MYSUBPLOT==1
            ylabel('Division location, S');
        end
    end
    xlim(  [0,max(allTimes)]  );
    MW_makeplotlookbetter(20);

end

%% Statistics on newborn cells
if any(strcmp(RUNSECTIONSFILADIV,{'all','newbornStats'}))

    allLengths = [myLengthNewborns{:}];
    
    hNBS = figure; clf; hold on;    
    
    [counts, edges] = histcounts(allLengths,100)
    dx=edges(2)-edges(1);
    centers = edges(2:end)-(dx)/2;

    pdf = counts./sum(counts).*dx;
    
    plot(centers, pdf,'LineWidth',3)
        
    
    MW_makeplotlookbetter(20);
    xlabel('Size of newborn cells (?m)');
    ylabel('Probability');


    %
    %===
    
    under5 = allLengths<5;
    [countsUnder5, edgesUnder5] = histcounts(allLengths(under5),200)
    dxU5 = edgesUnder5(2)-edgesUnder5(1);
    centersUnder5 = edgesUnder5(2:end)-(dxU5)/2;

    pdfUnder5 = countsUnder5./sum(countsUnder5).*dxU5;
    
    hNBSInset = figure; clf; hold on;% inset
    plot(centersUnder5,pdfUnder5,'LineWidth',3);
    
    set(gca,'YTick',[]); % note the counts are different w. main plot because the binsizes are different
    
    MW_makeplotlookbetter(20);
    xlim([1,3]);

    ymax = max(pdfUnder5)*1.1;
    ylim([0,ymax]);
    plot([1.78 1.78],[0,ymax],'LineWidth',3);
    
    disp(['Datapoints in distribution, N=' num2str(sum(counts))]);
    
end


%% Division preference for location?
if any(strcmp(RUNSECTIONSFILADIV,{'all','locationPreference'}))

    regimeNr = 3;
    regime = WINDOWBORDERS(regimeNr:regimeNr+1)
    allParentLengths = [myLengthParents{:}];
    allRatios = [Ratios{:}];

    selectedDivisionIndices = allParentLengths>regime(1) & allParentLengths<regime(2);
    selectedRatios = allRatios(selectedDivisionIndices);
    disp(['Based on N=' num2str(numel(selectedRatios)) ' division events.']);
    
    % figure
    hlocPref = figure; clf; hold on; [n,centers]=hist(selectedRatios,50); 
    
    % the bar plot
    %bar(centers,n,'FaceColor',linspecer(1),'EdgeColor','None');
    %barh(centers,n,'FaceColor',linspecer(1),'EdgeColor','None');
    barh(centers,n,'FaceColor',BARCOLORS(3,:),'EdgeColor','None');
    
    
    %xlim([0,0.5]);
    %myylim=max(n)*1.1; % for normal bars
    %ylim([0,myylim])   % for normal bars
    myxlim=max(n)*1.1; % horizontal bars
    xlim([0,myxlim]) % horizontal bars
    
    % boundary lines
    %plot([.33 .33],[0,myylim],':k','LineWidth',3); % for normal bars
    %plot([.66 .66],[0,myylim],':k','LineWidth',3); % for normal bars  
    plot([0,myxlim],[.33 .33],':k','LineWidth',2); % horizontal bars
    plot([0,myxlim],[.66 .66],':k','LineWidth',2);  % horizontal bars  

    % regime centers
    %plot([1/6,0.5,5/6],[0 0 0],'^k','MarkerFaceColor','k','MarkerSize',14); % for normal bars  
    plot([0 0 0],[1/6,0.5,5/6],'>k','MarkerFaceColor','k','MarkerSize',7); % horizontal bars

    % cosmetics
    %xlabel('Division location, S');
    %ylabel('Count');
    ylabel('Location S');
    xlabel('Count');
    yticks([0,1]);
    MW_makeplotlookbetter(20);
    
    % statistics
    site1=[0,0.33];
    site2=[0.33,0.66];
    site3=[0.66,1];
    indicesSite1 = selectedRatios<site1(2);
    indicesSite2 = selectedRatios>site2(1) & selectedRatios<site2(2);
    indicesSite3 = selectedRatios>site3(2);
    DivsInSite1 = selectedRatios(indicesSite1);
    DivsInSite2 = selectedRatios(indicesSite2);
    DivsInSite3 = selectedRatios(indicesSite3);
    n1=numel(DivsInSite1);
    n2=numel(DivsInSite2);
    n3=numel(DivsInSite3);

    % display statistics
    disp(['Stats this regime: ' num2str(n1) ' cells divide <.33 and ' num2str(n2) ' divide between .33-..66.']);
        
    

end

%% plotting derivative of length vs. cell size

if any(strcmp(RUNSECTIONSFILADIV,{'all','derivLvsL'}))

    %% plotting cell size against absolute growth rate

    allLenghts = []; alldVdt=[];
    for dataSetIdx=1:numel(datasetsPaths)

        schnitzcells = loadandrename(datasetsPaths{dataSetIdx});

        for schnitzIdx = 1:numel(schnitzcells)

            selectedFrames = schnitzcells(schnitzIdx).time>switchTimes(dataSetIdx);

            if numel(selectedFrames)>2

                %lengths = schnitzcells(schnitzIdx).length_fitNew(selectedFrames);
                %lengths = schnitzcells(schnitzIdx).length_skeleton(selectedFrames);
                lengths = schnitzcells(schnitzIdx).(LENGTHFIELD)(selectedFrames);
                t = schnitzcells(schnitzIdx).time(selectedFrames);

                dt = t(2:end)-t(1:end-1);
                dlengths = lengths(2:end)-lengths(1:end-1);

                dVdt = dlengths./dt;    
                lenghtAtPoints = lengths(1:end-1)+dlengths./2;

                allLenghts  = [allLenghts lenghtAtPoints];
                alldVdt     = [alldVdt dVdt];

            end

        end

        disp('Dataset done');
    end

    %%

    %NRBINS=20;
    NRBINS=5;
    %NRBINS=20;
    
    hdLdt1=figure(1); clf; hold on;
    scatter(allLenghts,alldVdt)

    edges=linspace(0,max(allLenghts),NRBINS);
    [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins]=binnedaveraging({allLenghts},{alldVdt},edges);

    errorbar(binCenters,meanValuesForBins,stdErrValuesForBins,'Color','k','LineWidth',4)
    plot(binCenters,meanValuesForBins,'-ok','LineWidth',4);

    ylim([-.1,.2]); xlim([0,40]);

    xlabel('L (µm)');
    ylabel(['dL/dt (µm/min)']);

    MW_makeplotlookbetter(20);

    %
    hdLdt2=figure(2); clf; hold on;
    errorbar(binCenters,meanValuesForBins,stdErrValuesForBins,'Color','k','LineWidth',2)
    ylim([-0,.1]); xlim([0,40]);

    xlabel('L (µm)');
    ylabel(['dL/dt (µm/min)']);

    MW_makeplotlookbetter(20);
    ylim([0,max(meanValuesForBins)*1.1]);
    
    disp(['N=' num2str(numel(allLenghts)) ' datapoints used']);

end
%% grandmother plot
%{
figure(11); clf; hold on;
for dataSetIndex = 1:numel(datasetsPaths)
        
    for i = 1:numel(myLengthSumNewborns{dataSetIndex})

        %currentLifeTime = lifeTimes{dataSetIndex}(myNewBornSchnitzNrs{dataSetIndex}(i)); % daughter lifetime
        currentTime = myLifeTimeParentsCorrectedQuickDiv{dataSetIndex}(i); % parent lifetime
        %currentTime = birthTimes{dataSetIndex}(i); % time at which born
        
        % determine color if lifetime known
        %{
        if isnan(currentTime)
            colorForThisDataPoint = [.5 .5 .5];
        else    
            colorForThisDataPoint = ...
                timeColormap(...
                    round(currentTime/fromLifeTimeToColorCodeFactor)+1 ...
                    ,:);
        end
        %}
    
        if myLifeTimeParentsCorrectedQuickDiv{dataSetIndex}(i)<20
            
            colorForThisDataPoint=[.7 .7 .7];
            theMarker = 'o';
            MarkerFaceColor='None';
            
        else
            
            if listWhichCorrectedQuickDiv{dataSetIndex}(i)==1
                colorForThisDataPoint='r';
                theMarker = 'o';
                MarkerFaceColor='r';
            else
                colorForThisDataPoint='b';
                theMarker = 'o';
                MarkerFaceColor='b';
            end                       
            
        end
        
        plot(myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}(i),RatiosCorrectedQuickDiv{dataSetIndex}(i),theMarker,'MarkerSize',10,...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',MarkerFaceColor);

        plot(myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}(i),1-RatiosCorrectedQuickDiv{dataSetIndex}(i),theMarker,'MarkerSize',10,...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',MarkerFaceColor);      
        
        
        % plot original data
        if listWhichCorrectedQuickDiv{dataSetIndex}(i)==1
            
            plot(myLengthSumNewborns{dataSetIndex}(i),Ratios{dataSetIndex}(i),'s','MarkerSize',7,...
                        'MarkerEdgeColor',[1 .5 0],...
                        'MarkerFaceColor','None',...
                        'LineWidth',1);

            plot(myLengthSumNewborns{dataSetIndex}(i),1-Ratios{dataSetIndex}(i),'s','MarkerSize',7,...
                        'MarkerEdgeColor',[1 .5 0],...
                        'MarkerFaceColor','None',...
                        'LineWidth',1);
                    
        end
                
    end
    
end

xlabel('Summed daughter length [um]');
ylabel('L_d/L_p');


%{
hc = colorbar;
%ax = gca;
%ax.XTickLabel = {'1000','1'};
caxis([0,longestTime])
ylabel(hc, 'Daughter interdivision time')
%}

MW_makeplotlookbetter(20);

xlim([0,40]);

disp('section done');
%}
%% Plot historgrams per window again, now with corrected/uncorrected dataset
%{
NRREGIMES = 4;
WINDOWBORDERS = [5,10,17,20,30];
WINDOWBORDERS = [2,9,16,23,30];
WINDOWBORDERS = [2:7:30];
WINDOWBORDERS = [3:6:30];
LINEWIDTH =2;

% figure stuff
h=figure(12); clf; hold on
myColors = linspecer(numel(WINDOWBORDERS)-1);

% calculate params
mybins = linspace(0,1,HISTNRBINS);
mycenters = mybins(2:end)-(mybins(2:end)-mybins(1:end-1))/2;

for dataSetIndex = 1:numel(datasetsPaths)
    %dataSetIndex=1; % TEMP REMOVE   

    for windowIndex = 1:(numel(WINDOWBORDERS)-1)

        windowLeft = WINDOWBORDERS(windowIndex);
        windowRight  = WINDOWBORDERS(windowIndex+1);
        windowSize = windowRight-windowLeft;

        % get data
        currentWindowIndices = ...
            find((myLengthSumNewborns{dataSetIndex}>windowLeft) == (myLengthSumNewborns{dataSetIndex}<windowRight));
        currentTotalLength = myLengthSumNewborns{dataSetIndex}(currentWindowIndices);
        currentRatios = Ratios{dataSetIndex}(currentWindowIndices);
        
        currentWindowIndicesCorrectedQuickDiv = ...
            find((myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}>windowLeft) == (myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}<windowRight));        
        currentTotalLengthCorrectedQuickDiv = myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}(currentWindowIndicesCorrectedQuickDiv);
        currentRatiosCorrectedQuickDiv = RatiosCorrectedQuickDiv{dataSetIndex}(currentWindowIndicesCorrectedQuickDiv);

        % top row
        % ===
        subplot(2,1,1); hold on
        plot(currentTotalLength,currentRatios,'o','Color',myColors(windowIndex,:),'MarkerFaceColor','None')
        plot(currentTotalLengthCorrectedQuickDiv,currentRatiosCorrectedQuickDiv,'o','Color',myColors(windowIndex,:),'MarkerFaceColor',myColors(windowIndex,:))
        % cosmetics
        ylim([0,1]);
        xlim([LEFTX,RIGHTX]);
        xlabel(['Parent cell size by ' LENGTHFIELD],'Interpreter','None');
        ylabel('L_{child}/L_{parent}');
        MW_makeplotlookbetter(15)
        title([WHATDATA ' condition']);
        
        % bottom row with hists
        % ===
        subplot(2,1,2); hold on;    
        % predicted location        
        nrLocations = windowIndex;
        for j = 1:2:(nrLocations*2-1)
            plot([windowLeft, windowRight], [(j)/(2*nrLocations) (j)/(2*nrLocations)],'-','Color',[.5 .5 .5],'LineWidth',2);
        end        
        % histograms
        [counts,centers]=hist(currentRatios,mycenters);        
        modcounts = (windowSize*counts/sum(counts)+windowLeft);
        plot(modcounts,centers,':','Color',myColors(windowIndex,:),'LineWidth',1)
        
        [countsCorrectedQuickDiv,centers]=hist(currentRatiosCorrectedQuickDiv,mycenters);
        modcountsCorrectedQuickDiv = (windowSize*countsCorrectedQuickDiv/sum(countsCorrectedQuickDiv)+windowLeft);
        plot(modcountsCorrectedQuickDiv,centers,'-','Color',myColors(windowIndex,:),'LineWidth',LINEWIDTH)
        
        % cosmetics
        ylim([0,1]);
        xlim([LEFTX,RIGHTX]);
        xlabel(['Histograms (normalized)']);
        ylabel('L_{child}/L_{parent}');
        MW_makeplotlookbetter(15)
        set(gca,'xtick',[])

    end
end

if ~NOSAVEPLEASE
    saveas(h, [PLOTSAVEDIR WHATDATA '_ratios_histograms.fig']);
    saveas(h, [PLOTSAVEDIR WHATDATA '_ratios_histograms.tif']);
end
%}
%%
%{
figure(FIGURENUMBERS(end)+3); clf; hold on;

plot(histArea(2,:),histArea(1,:),'-r');
plot(histSkel(2,:),histSkel(1,:),'-b');
%}



%% 
if exist('ONESTOTAKE','var')
    disp(['Whole condition ' WHATDATA ' done (#' num2str(ONESTOTAKE) ').']);
else
    disp(['Whole condition ' WHATDATA ' done.']);
end





