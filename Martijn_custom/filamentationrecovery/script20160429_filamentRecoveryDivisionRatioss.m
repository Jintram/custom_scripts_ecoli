
%% 

% NOTE ABOUT CorrectedQuickDiv
% ===
% Note that I tried to correct for cells that divide fast, but that this is
% impossible, since we don't know the orientation of the cells.

% Plotting division ratios
%WHATDATA = 'sulA';
WHATDATA = 'temperature';
%WHATDATA = 'simulated';

USESYMMETRY=1;
HISTNRBINS=50;

FIGURENUMBERS=[1 1 1 1 1 1 1]; % TODO this option does not work any more

%PLOTCOLORS = [0 0 1;0 0 1;0 0 1];
%PLOTCOLORS = [0,0,139 ; 0,191,255  ; 	135,206,250]./255; % shades of blue
%PLOTCOLORS = [255, 0, 0; 204, 0, 204 ; 255, 102, 0]./255; % shades of blue
PLOTCOLORS = linspecer(7);

%LENGTHFIELD = 'areaPixels';
%LENGTHFIELD = 'length_fitNew';
LENGTHFIELD = 'length_skeleton';
%LENGTHFIELD = 'cellLengths';

TIMEFIELD = 'time';
%TIMEFIELD = 'times';

if ~exist('PLOTSAVEDIR','var') & ~exist('NOSAVEPLEASE','var');
    error('Set PLOTSAVEDIR please. Or set NOSAVEPLEASE');
end
if ~exist('NOSAVEPLEASE','var')
    NOSAVEPLEASE=0;
end

%{
% Please execute this code to define plot dir and make log of scriptname
PLOTSAVEDIR = 'U:\PROJECTS\B_filamentationRecovery\Contributions_Martijn\figures_fig_files\';
makenoteinreadme(PLOTSAVEDIR,'script20160429_filamentRecoveryDivisionRatioss');
%}

%% The sulA datasets

if strcmp(WHATDATA, 'sulA')
    datasetsPaths = ...
        { ...
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos1crop\data\pos1crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos2crop\data\pos2crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos3crop\data\pos3crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos4crop\data\pos4crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos7crop\data\pos7crop-Schnitz.mat'...
        }
elseif strcmp(WHATDATA, 'temperature')
    datasetsPaths = ...
        { ...
        'G:\EXPERIMENTAL_DATA_2016\c_completely_analyzed\2016-03-23\pos4crop\data\pos4crop-Schnitz.mat',...
        ['G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos2crop\data\pos2crop-Schnitz.mat'],...
        }
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
%% 
for i = FIGURENUMBERS
    figure(i); clf; hold on;
end

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

myLifeTimeParents = {}; 
myLifeTimesSchnitzes = {};
allLengths = {};
%myLengthParentsCorrectedQuickDiv={};
%myLifeTimeParentsCorrectedQuickDiv={};
%myLengthSumNewbornsCorrectedQuickDiv={};
%listWhichCorrectedQuickDiv={};

for dataSetIndex = 1:numel(datasetsPaths)
    
    if ~exist('simulatedschnitzcells','var')
        load(datasetsPaths{dataSetIndex});
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
    
    for i = 1:numel(schnitzcells)

        LengthParent = schnitzcells(i).(LENGTHFIELD)(end);
        
        daughterSchnitz1 = schnitzcells(i).D;
        daughterSchnitz2 = schnitzcells(i).E;

        if ~any([daughterSchnitz1,daughterSchnitz2]==0)

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
                myLifeTimeParents{dataSetIndex}(end+1) = schnitzcells(i).interDivTime; % daughter 1
                myLifeTimeParents{dataSetIndex}(end+1) = schnitzcells(i).interDivTime; % daughter 2
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
    if isfield(schnitzcells,'interDivTime')
        myLifeTimesSchnitzes{dataSetIndex} = [schnitzcells.interDivTime];
    end
    for i = 1:numel(schnitzcells)
        thisSchnitzLengths = schnitzcells(i).(LENGTHFIELD);
        allLengths{dataSetIndex}(i) = thisSchnitzLengths(1);
        birthTimes{dataSetIndex}(i) = schnitzcells(i).(TIMEFIELD)(1);
    end
    
    if ~isfield(schnitzcells,'interDivTime')
        warning('interDivTime field was not found')
    end
    
end

disp('Section done');

%% Now calculate ratios and plot them
if ~exist('SANITYCHECKSYMMETRY','var')
    SANITYCHECKSYMMETRY=1;
end

Ratios={}; 
% RatiosCorrectedQuickDiv = {};
for dataSetIndex = 1:numel(datasetsPaths)
    
    %%
    % user given parameters    
    LEFTX =  0;   

    % calculated parameters
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

    % create figure
    figure(1); 

    % plot helping lines at 1/2n
    % for i=1:5
    %     plot([rightx, LEFTX], [.5/i .5/i],':','Color',[.5 .5 .5],'LineWidth',2)
    %     plot([rightx, LEFTX], 1-[.5/i .5/i],':','Color',[.5 .5 .5],'LineWidth',2)
    % end
    N=5;
    for i=1:N
        for j = 1:(i*2-1)
            plot([0, RIGHTX], [(j)/(2*i) (j)/(2*i)],'-','Color',[.5 .5 .5],'LineWidth',N-i+1)
        end
    end

    % target length line
    x = (LONGESTNORMALDIVSIZEPARENT):RIGHTX;
    y = (LONGESTNORMALDIVSIZEPARENT/2)./x;
    plot(x,y,'--','Color','k','LineWidth',2);    
    
    % plot ratios
    Ratios{dataSetIndex} = myLengthNewborns{dataSetIndex}./myLengthSumNewborns{dataSetIndex};
    %RatiosCorrectedQuickDiv{dataSetIndex} = myLengthNewborns{dataSetIndex}./myLengthSumNewbornsCorrectedQuickDiv{dataSetIndex};
    plot(myLengthSumNewborns{dataSetIndex},Ratios{dataSetIndex},'o', 'Color', PLOTCOLORS(dataSetIndex,:),'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',PLOTCOLORS(dataSetIndex,:));
    if SANITYCHECKSYMMETRY    
        plot(myLengthSumNewborns{dataSetIndex},1-Ratios{dataSetIndex},'x', 'Color', [.5 .5 .5],'LineWidth',2,'MarkerSize',15);
    end

    % cosmetics
    ylim([0,1]);
    xlim([LEFTX,RIGHTX]);

    xlabel(['Parent cell size by ' LENGTHFIELD],'Interpreter','None');
    ylabel('L_{child}/L_{parent}');

    MW_makeplotlookbetter(15)

end

disp('section done');

%% Now calculate overall histogram
h=figure(2); clf; hold on;

% calculate hist
[count,bincenters] = hist([Ratios{:}],HISTNRBINS);

% helping lines
N=5; highestcount=max(count);
for i=1:N
    for j = 1:(i*2-1)
        plot([(j)/(2*i) (j)/(2*i)],[0, highestcount],'-','Color',[.5 .5 .5],'LineWidth',N-i+1)
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

%% Sanity check
% load(datasetsPaths{2})
figure(3); clf; hold on

% x=y line
plot([1,10^6],[1,10^6],'-k')

% data
for i=1:numel(datasetsPaths)
    plot(myLengthSumNewborns{dataSetIndex}, myLengthParents{dataSetIndex},'x')
end

% cosmetics
axis equal;
xlim([0 max([myLengthSumNewborns{:}])]);
ylim([0 max([myLengthParents{:}])]);
xlabel('Summed length newborns');
ylabel('Length of parent');

%% Plot historgrams per window
NRREGIMES = 4;
%WINDOWBORDERS = [5,10,17,20,30];
%WINDOWBORDERS = [2,9,16,23,30];
%WINDOWBORDERS = [2:7:30];
WINDOWBORDERS = [3:6:30];
regionsDouble = {[0,0.5],[0,0.5],[0,2/6,.5],[0,2/8,.5]}; % boundaries around matching ratios
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


if ~NOSAVEPLEASE
    saveas(h, [PLOTSAVEDIR WHATDATA '_ratios_histograms.fig']);
    saveas(h, [PLOTSAVEDIR WHATDATA '_ratios_histograms.tif']);
end

%% Plot simpler figure w. ratios
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

%% Plot lifetime against birth length

h=figure(5); clf; hold on
myColors = linspecer(numel(datasetsPaths));
myPlotMarkers = 'os^vd<>';

title([WHATDATA ' condition']);

%
meanInterDivisionTimes=[];
for dataSetIndex = 1:numel(datasetsPaths)

    plot(allLengths{dataSetIndex},myLifeTimesSchnitzes{dataSetIndex},'.',...
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

    end

    % plot
    %errorbar(meansLengths,meansLifeTimes,stdLifeTimes,'ok','MarkerFaceColor','k')
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

%% Plot with time-coded ratios
MINUTELIFETIMETRESHOLD=20;

h1=figure(7); clf; hold on
h2=figure(8); clf; hold on

% Get a colormap (is 64 long by default)
timeColormap = colormap(jet); % winter
rowsInColormap = size(timeColormap,1);

% 
longestTime = max([myLifeTimesSchnitzes{:}]);
%longestTime = max([birthTimes{:}]);
fromLifeTimeToColorCodeFactor = longestTime/(rowsInColormap-1);

for dataSetIndex = 1:numel(datasetsPaths)

    for i = 1:numel(myLengthSumNewborns{dataSetIndex})

        %currentLifeTime = lifeTimes{dataSetIndex}(myNewBornSchnitzNrs{dataSetIndex}(i)); % daughter lifetime
        currentTime = myLifeTimeParents{dataSetIndex}(i); % parent lifetime
        %currentTime = birthTimes{dataSetIndex}(i); % time at which born
        
        % determine color if lifetime known
        if isnan(currentTime)
            colorForThisDataPoint = [.5 .5 .5];
        else    
            colorForThisDataPoint = ...
                timeColormap(...
                    round(currentTime/fromLifeTimeToColorCodeFactor)+1 ...
                    ,:);
        end

        % plot all datapoint to one figure
        figure(h1.Number);
        plot(myLengthSumNewborns{dataSetIndex}(i),Ratios{dataSetIndex}(i),'o','MarkerSize',10,...
                'MarkerEdgeColor',colorForThisDataPoint,...
                'MarkerFaceColor',colorForThisDataPoint)
        if USESYMMETRY
            plot(myLengthSumNewborns{dataSetIndex}(i),1-Ratios{dataSetIndex}(i),'o','MarkerSize',10,...
                'MarkerEdgeColor',colorForThisDataPoint,...
                'MarkerFaceColor',colorForThisDataPoint)
        end
        
        
        % plot selection of datapoints to 2nd figure
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
    end
    
end

for myfignum = h1.Number:h2.Number
    figure(myfignum);

    xlabel('Summed daughter length [um]');
    ylabel('L_d/L_p');
    
    hc = colorbar;
    %ax = gca;
    %ax.XTickLabel = {'1000','1'};
    caxis([0,longestTime])
    ylabel(hc, 'Daughter interdivision time')

    MW_makeplotlookbetter(20);
end

if exist('MYXLIM','var')
    figure(7); xlim(MYXLIM); figure(8); xlim(MYXLIM);
else
    disp('Set MYXLIM to [limleft limright] to adjust plots x axis to custom values');
end

if ~NOSAVEPLEASE
    saveas(h1, [PLOTSAVEDIR WHATDATA '_ratiosByInterdiv.fig']);
    saveas(h1, [PLOTSAVEDIR WHATDATA '_ratiosByInterdiv.tif']);

    saveas(h2, [PLOTSAVEDIR WHATDATA '_ratiosByInterdivSelected.fig']);
    saveas(h2, [PLOTSAVEDIR WHATDATA '_ratiosByInterdivSelected.tif']);
end

%% Histogram of lifetimes
h=figure(9); clf; hold on

hist([myLifeTimeParents{:}],30)
xlabel('Interdivision lifetime [min]');
ylabel('Count')

MW_makeplotlookbetter(20);

if ~NOSAVEPLEASE
    saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_generalhistogram.fig']);
    saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_generalhistogram.tif']);
end

%% histogram of cell sizes

figure(10); clf; hold on;
allLengths = [schnitzcells(:).(LENGTHFIELD)];
histogram(allLengths,100);

xlabel('Length (um)');
ylabel('Count');
title('Cell lenghts observed during movie');
MW_makeplotlookbetter(15);

%% Plot with time-coded ratios
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

    for i = 1:numel(myLengthSumNewborns{dataSetIndex})

        %currentLifeTime = lifeTimes{dataSetIndex}(myNewBornSchnitzNrs{dataSetIndex}(i)); % daughter lifetime
        currentTime = myLifeTimeParents{dataSetIndex}(i); % parent lifetime
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
        plot(myLengthSumNewborns{dataSetIndex}(i),Ratios{dataSetIndex}(i),'o','MarkerSize',10,...
                'MarkerEdgeColor',colorForThisDataPoint,...
                'MarkerFaceColor',faceForPoint)
        if USESYMMETRY
            plot(myLengthSumNewborns{dataSetIndex}(i),1-Ratios{dataSetIndex}(i),'o','MarkerSize',10,...
                'MarkerEdgeColor',colorForThisDataPoint,...
                'MarkerFaceColor',faceForPoint)
        end
                
    end
    
end

xlabel('Summed daughter length [um]');
ylabel('Relative division location (S)');

MW_makeplotlookbetter(20);

if exist('MYXLIM','var')
    figure(h1.Number); xlim(MYXLIM);
else
    disp('Set MYXLIM to [limleft limright] to adjust plots x axis to custom values');
end

if ~NOSAVEPLEASE
    saveas(h1, [PLOTSAVEDIR WHATDATA '_newRutgerPlot.fig']);
    saveas(h1, [PLOTSAVEDIR WHATDATA '_newRutgerPlot.tif']);
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
disp('Done.');




