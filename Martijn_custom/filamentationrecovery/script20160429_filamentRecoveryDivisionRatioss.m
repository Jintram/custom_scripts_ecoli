
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
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos1crop\data\pos1crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos2crop\data\pos2crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos3crop\data\pos3crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos4crop\data\pos4crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos7crop\data\pos7crop-Schnitz.mat'...
        }
    switchTimes = [0 0 0 0 0];
    SPECIALCASE=0;
elseif strcmp(WHATDATA, 'temperature')
    datasetsPaths = ...
        { ...
        'G:\EXPERIMENTAL_DATA_2016\c_completely_analyzed\2016-03-23\pos4crop\data\pos4crop-Schnitz.mat',...
            ... No fluor signal
        ['G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos2crop\data\pos2crop-Schnitz.mat'],...
            ... 
        }
    switchTimes = [450, 329];
    SPECIALCASE=0;
elseif strcmp(WHATDATA, 'tetracycline')
    
    % Load data from gathered file:
    %SPECIALCASE=1;
    %load('G:\FilamentationRecoveryData\Dmitry\3_Rutgers_project\Recent Data\rutgers old data\DE analysis\analysis of clusters\all_rutgers_divisions_data_gen_shifted.mat');    
    %datasetsPaths=exp_names_all;
    
    % Or from Rutger's files
    ONESTOTAKE=[1:5]; % ONESTOTAKE=[1:5]; ONESTOTAKE=[6:8]; ONESTOTAKE=[9:11];
    disp(['Note that certain datasets are selected (' num2str(ONESTOTAKE) ').']);
    switchTimes = ...
        [1373.5417      587.66667      261.13333        741.325      252.38333      987.86667        1017.95         997.05      881.90833      1018.7083         1394.6];
        % determined by 
        % script20161222_branchplottingrutgerdata
        % which calculates time at which growth rate is restored halfway
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

    myLifeTimeParents = {}; 
    myLifeTimesSchnitzes = {};
    allLengths = {};
    birthTimes = {};
    %myLengthParentsCorrectedQuickDiv={};
    %myLifeTimeParentsCorrectedQuickDiv={};
    %myLengthSumNewbornsCorrectedQuickDiv={};
    %listWhichCorrectedQuickDiv={};
    myCellLengthsPerSchnitz = {}; 
    myAddedLengthPerSchnitz = {};

    for dataSetIndex = 1:numel(datasetsPaths)

        if ~exist('simulatedschnitzcells','var')
            if ~SPECIALCASE
                schnitzcells = loadandrename(datasetsPaths{dataSetIndex});            
            else
                schnitzcells=S_all_shifted{dataSetIndex};
            end
        end

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
        %myLengthParentsCorrectedQuickDiv{dataSetIndex}=[];
        %myLifeTimeParentsCorrectedQuickDiv{dataSetIndex}=[];
        %myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex}=[];
        %listWhichCorrectedQuickDiv{dataSetIndex}=[];
        
        for binIdx = 1:numel(schnitzcells)

            LengthParent = schnitzcells(binIdx).(LENGTHFIELD)(end);

            daughterSchnitz1 = schnitzcells(binIdx).D;
            daughterSchnitz2 = schnitzcells(binIdx).E;            
            
            if ~any([daughterSchnitz1,daughterSchnitz2]==0)
                
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
    hS=scatter(selectedXdata,selectedYdata,'filled',...
            'MarkerFaceColor',markerColor,'MarkerEdgeColor','none','MarkerFaceAlpha',1);
    set(hS, 'SizeData', markerWidth^2);

    errorbar(binCentersDynData,meanValuesForBinsDynData,stdValuesForBinsDynData,'ok','LineWidth',3,'MarkerFaceColor','k');

    MW_makeplotlookbetter(20);
    xlabel('Birth size (\mum)');
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
            legendTitles{end+1} = ['' sprintf('%0.1f', binCentersDynData(binIdx)) '\mum (N=' num2str(totalCount) ')'];
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
    
    %%
    hFig2D=figure(103); clf; hold on;

    markerColor=[.5 .5 .5];

    alldatax=[];alldatay=[];
    for datasetIdx=1:numel(datasetsPaths)

        bornAfterSwitchIdxs = (birthTimes{datasetIdx}-switchTimes(datasetIdx))>0;    
        selectedLifeTimeData = myLifeTimesSchnitzes{datasetIdx}(bornAfterSwitchIdxs);

        baseColors = linspecer(3);
        myColors = makeColorMap(baseColors(1,:), baseColors(2,:), baseColors(3,:));
        lifeTimesNormalized=round(...
            (selectedLifeTimeData-min(selectedLifeTimeData))./    ...
            (max(selectedLifeTimeData)-min(selectedLifeTimeData)) ...
            .*99 ...
            )+1;
        lifeTimesNormalized(isnan(lifeTimesNormalized))=1;
        markerColors=myColors(lifeTimesNormalized,:);

        l=scatter(  myCellLengthsPerSchnitz{datasetIdx}(bornAfterSwitchIdxs),...
                    myAddedLengthPerSchnitz{datasetIdx}(bornAfterSwitchIdxs),...
                    [],markerColors,'filled','MarkerFaceAlpha',.5)

        %l=scatter(  myCellLengthsPerSchnitz{datasetIdx}(bornAfterSwitchIdxs),...
        %    myAddedLengthPerSchnitz{datasetIdx}(bornAfterSwitchIdxs));
        %set(l,'MarkerFaceColor',markerColors,'MarkerEdgeColor','none','MarkerFaceAlpha',.5);

        alldatax = [alldatax myCellLengthsPerSchnitz{datasetIdx}(bornAfterSwitchIdxs)];
        alldatay = [alldatay myAddedLengthPerSchnitz{datasetIdx}(bornAfterSwitchIdxs)];
        
    end

    [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=binnedaveraging({alldatax},{alldatay},linspace(0,20,10));
    errorbar(binCenters,meanValuesForBins,stdValuesForBins,'ok','LineWidth',2);
    
    %
    ylim([0,4]); 
    xlim([0,20]);
    xlabel('Cell length'); ylabel('Added length \Delta (\mum)');
    MW_makeplotlookbetter(20);

    colormap(myColors);
    colorbar;

    %% Now also create pdfs
    xEdges=linspace(0,20,10);
    xCenters = xEdges(1:end-1)+(xEdges(2)-xEdges(1))/2;
    yEdges=linspace(0,4,10);
    dY = (yEdges(2)-yEdges(1));
    yCenters = yEdges(1:end-1)+dY/2;
    [pdfsForBins, binCenters]=binnedpdfs({alldatax},{alldatay},xEdges,yEdges);
    
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
        
        legendTitles{end+1} = ['L_b=' sprintf('%0.1f', xCenters(binIdx)) '\mum (N=' num2str(totalN) ')'];
        
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

    for dataSetIndex = 1:numel(datasetsPaths)

        for binIdx = 1:numel(myLengthSumNewborns{dataSetIndex})

            %currentLifeTime = lifeTimes{dataSetIndex}(myNewBornSchnitzNrs{dataSetIndex}(i)); % daughter lifetime
            currentTime = myLifeTimeParents{dataSetIndex}(binIdx); % parent lifetime
            %currentTime = birthTimes{dataSetIndex}(i); % time at which born

            % determine color if lifetime known
            if isnan(currentTime) || ~COLORCODEMARKERS
                colorForThisDataPoint = [.7 .7 .7];
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
                    'MarkerFaceAlpha',0.15,...
                    'MarkerEdgeAlpha',0.15);
            if USESYMMETRY
                scatter(myLengthSumNewborns{dataSetIndex}(binIdx),1-Ratios{dataSetIndex}(binIdx),12^2,'filled',...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',colorForThisDataPoint,...
                    'MarkerFaceAlpha',0.15,...
                    'MarkerEdgeAlpha',0.15);
            end


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

        xlabel('Length of mother cell [\mum]');
        ylabel('Relative division location, S');

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
        N=5;
        for windowIndex=1:(numel(WINDOWBORDERS)-1)
            for j = 1:2:(windowIndex*2-1)
                plot([WINDOWBORDERS(windowIndex), WINDOWBORDERS(windowIndex+1)], [(j)/(2*windowIndex) (j)/(2*windowIndex)],'-','Color',BARCOLORS(windowIndex,:),'LineWidth',8)
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

        for binIdx = 1:numel(myLengthSumNewborns{dataSetIndex})

            %currentLifeTime = lifeTimes{dataSetIndex}(myNewBornSchnitzNrs{dataSetIndex}(i)); % daughter lifetime
            currentTime = myLifeTimeParents{dataSetIndex}(binIdx); % parent lifetime
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
            plot(myLengthSumNewborns{dataSetIndex}(binIdx),Ratios{dataSetIndex}(binIdx),'o','MarkerSize',10,...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',faceForPoint)
            if USESYMMETRY
                plot(myLengthSumNewborns{dataSetIndex}(binIdx),1-Ratios{dataSetIndex}(binIdx),'o','MarkerSize',10,...
                    'MarkerEdgeColor',colorForThisDataPoint,...
                    'MarkerFaceColor',faceForPoint)
            end

        end

    end

    xlabel('Summed daughter length [um]');
    ylabel('Relative division location (S)');

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

end
%% Size time evolution plot
if any(strcmp(RUNSECTIONSFILADIV,{'all','sizeTraces'}))

    hSizeTime = figure(12); clf; hold on;

    MARKERSIZE=8;
    LINEWIDTH=6;
    FONTSIZE=20;
    SUBDIR = 'trace_dynamics\';
    NRSTARTSCHNITZES = 1;
    COLORDOTS = 0; % 1 or 0
    MAXTRACESHIGHLIGHTED=4;

    traceColors = linspecer(NRSTARTSCHNITZES*numel(datasetsPaths));

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
        markerWidth = 8;
        hS=scatter(timeScatterData-switchTimes(dataSetIndex),lengthScatterData,'filled',...
                'MarkerFaceColor',markerColor,'MarkerEdgeColor','none','MarkerFaceAlpha',1);
        set(hS, 'SizeData', markerWidth^2);

        % Now get ends of schnitzcells lineages
        endSchnitzesIdxs = find([schnitzcells(:).D]==0);
        theEndTimes = arrayfun(@(x) schnitzcells(x).(TIMEFIELD)(end), endSchnitzesIdxs);
        theEndLengths = arrayfun(@(x) schnitzcells(x).(LENGTHFIELD)(end), endSchnitzesIdxs);
        endOfLineages{dataSetIndex}= [theEndTimes; theEndLengths];

        % Now also make example traces
        [myCurrentTraces, myCurrentDivisions] = MW_gettracefromschnitzcells(schnitzcells,NRSTARTSCHNITZES,TIMEFIELD,LENGTHFIELD,'longest');
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

    [myShortTrace, myShortTraceDivisions] = MW_gettracefromschnitzcells(schnitzcells,NRSTARTSCHNITZES,TIMEFIELD,LENGTHFIELD,'shortest');
    % now plot traces that were made
    % following short schnitzes
    startSchnitzIdx=1;
    plot(myShortTrace{startSchnitzIdx}(:,1)-switchTimes(dataSetIndex),myShortTrace{startSchnitzIdx}(:,2),'Color','k','LineWidth',LINEWIDTH);
    plot(myShortTraceDivisions{startSchnitzIdx}(:,1)-switchTimes(dataSetIndex),myShortTraceDivisions{startSchnitzIdx}(:,2),'o','Color','k','MarkerSize',MARKERSIZE,'MarkerFaceColor','k');
    % following long schnitzes
    counter=0;
    for dataSetIndex = 1:min(numel(datasetsPaths),MAXTRACESHIGHLIGHTED)
    for startSchnitzIdx = 1:NRSTARTSCHNITZES
        counter=counter+1;
        plot(myTraces{dataSetIndex}{startSchnitzIdx}(:,1)-switchTimes(dataSetIndex),myTraces{dataSetIndex}{startSchnitzIdx}(:,2),'Color',traceColors(counter,:),'LineWidth',LINEWIDTH);
        plot(myDivisions{dataSetIndex}{startSchnitzIdx}(:,1)-switchTimes(dataSetIndex),myDivisions{dataSetIndex}{startSchnitzIdx}(:,2),'o','Color',traceColors(counter,:),'MarkerSize',MARKERSIZE,'MarkerFaceColor',traceColors(counter,:));
    end
    end

    % cosmetics
    MW_makeplotlookbetter(FONTSIZE);
    xlabel('Time (min)');
    ylabel('Cell length (um)');
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
    for startSchnitzIdx = 1:NRSTARTSCHNITZES    
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





