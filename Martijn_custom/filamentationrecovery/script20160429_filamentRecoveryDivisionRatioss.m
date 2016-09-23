
%% 
% Plotting division ratios
WHATDATA = 'sulA';
WHATDATA = 'temperature';
%WHATDATA = 'simulated';

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
    % now there already should be a simulated schnitzcells
    datasetsPaths = ...
        { ...
        ...'D:\Local_Software\Martijn_extensions\Martijn_custom\filamentationrecovery\simulatedSchnitzcells\schnitz1_naivemodel_fastrecharge.mat', ...
        'D:\Local_Software\Martijn_extensions\Martijn_custom\filamentationrecovery\simulatedSchnitzcells\schnitz1_naivemodel_longrecharge.mat' ...
        };
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
% - myLengthNewborns{datasetIdx}(schnitz) 
% - myLengthParents{datasetIdx}(schnitz)
% - myLengthSumNewborns{datasetIdx}(schnitz)
%
% Where myLengthNewborns gives the length of a newborn, and myLengthParents
% gives the length of the corresponding parent. Because there are some
% subtleties (growth, geometrical distortion poles), myLengthSumNewborns is
% a better measure of the total length before division, and is simply the
% summed length of the two daughter cells.

lifeTimes = {}; 
allLengths = {};

for datasetIdx = 1:numel(datasetsPaths)
    
    load(datasetsPaths{datasetIdx});

    %% 
    figureIndex=FIGURENUMBERS(datasetIdx);

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
    myLengthNewborns{datasetIdx} = []; myLengthParents{datasetIdx} = []; myLengthSumNewborns{datasetIdx}=[];
    myNewBornSchnitzNrs{datasetIdx} = []; %myDaughter1schnitzNrs{datasetIdx} = []; myDaughter2schnitzNrs{datasetIdx} = [];
    
    for i = 1:numel(schnitzcells)

        LengthParent = schnitzcells(i).(LENGTHFIELD)(end);
        
        daughterSchnitz1 = schnitzcells(i).D;
        daughterSchnitz2 = schnitzcells(i).E;

        if ~any([daughterSchnitz1,daughterSchnitz2]==0)

            lengthDaughterSchnitz1 = schnitzcells(daughterSchnitz1).(LENGTHFIELD)(1);
            lengthDaughterSchnitz2 = schnitzcells(daughterSchnitz2).(LENGTHFIELD)(1);            
            
            % daughter 1
            myLengthNewborns{datasetIdx}(end+1) =     lengthDaughterSchnitz1;
            myLengthParents{datasetIdx}(end+1) =      LengthParent;
            myLengthSumNewborns{datasetIdx}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2;
            myNewBornSchnitzNrs{datasetIdx}(end+1) =  daughterSchnitz1;
            % daughter 2
            myLengthNewborns{datasetIdx}(end+1) =     lengthDaughterSchnitz2;
            myLengthParents{datasetIdx}(end+1) =      LengthParent;
            myLengthSumNewborns{datasetIdx}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2;
            myNewBornSchnitzNrs{datasetIdx}(end+1) =  daughterSchnitz2;
            
            % For later reference, make lookup table 
            %myDaughter1schnitzNrs{datasetIdx}(end+1) = daughterSchnitz1;
            %myDaughter2schnitzNrs{datasetIdx}(end+1) = daughterSchnitz2;                        
        end

    end

    %% Gather data on interdivision time vs. lifetime of cell
    if isfield(schnitzcells,'interDivTime')
        lifeTimes{datasetIdx} = [schnitzcells.interDivTime];    
    end
    birthTimes{datasetIdx} = NaN(1,numel(schnitzcells));
    allLengths{datasetIdx} = NaN(1,numel(schnitzcells));
    for i = 1:numel(schnitzcells)
        thisSchnitzLengths = schnitzcells(i).(LENGTHFIELD);
        allLengths{datasetIdx}(i) = thisSchnitzLengths(1);
        birthTimes{datasetIdx}(i) = schnitzcells(i).(TIMEFIELD)(1);
    end
    
end

%% Now calculate ratios and plot them
SANITYCHECKSYMMETRY=1;

Ratios={};
for datasetIdx = 1:numel(datasetsPaths)
    
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
        RIGHTX  = max(myLengthParents{datasetIdx})*1.1;
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
    Ratios{datasetIdx} = myLengthNewborns{datasetIdx}./myLengthSumNewborns{datasetIdx};
    plot(myLengthSumNewborns{datasetIdx},Ratios{datasetIdx},'x', 'Color', PLOTCOLORS(datasetIdx,:),'LineWidth',2);
    if SANITYCHECKSYMMETRY    
        plot(myLengthSumNewborns{datasetIdx},1-Ratios{datasetIdx},'o', 'Color', PLOTCOLORS(datasetIdx,:),'LineWidth',2);
    end

    % cosmetics
    ylim([0,1]);
    xlim([LEFTX,RIGHTX]);

    xlabel(['Parent cell size by ' LENGTHFIELD],'Interpreter','None');
    ylabel('L_{child}/L_{parent}');

    MW_makeplotlookbetter(15)

end

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
    plot(myLengthSumNewborns{datasetIdx}, myLengthParents{datasetIdx},'x')
end

% cosmetics
axis equal;
xlim([0 max([myLengthSumNewborns{:}])]);
ylim([0 max([myLengthParents{:}])]);
xlabel('Summed length newborns');
ylabel('Length of parent');

%% Plot historgrams per window
NRREGIMES = 4;
WINDOWBORDERS = [5,10,17,20,30];
WINDOWBORDERS = [2,9,16,23,30];
WINDOWBORDERS = [2:7:30];
WINDOWBORDERS = [3:6:30];
LINEWIDTH =2;

% figure stuff
h=figure(4); clf; hold on
myColors = linspecer(numel(WINDOWBORDERS)-1);

% calculate params
mybins = linspace(0,1,HISTNRBINS);
mycenters = mybins(2:end)-(mybins(2:end)-mybins(1:end-1))/2;

for datasetIdx = 1:numel(datasetsPaths)
    %datasetIdx=1; % TEMP REMOVE   

    for windowIndex = 1:(numel(WINDOWBORDERS)-1)

        windowLeft = WINDOWBORDERS(windowIndex);
        windowRight  = WINDOWBORDERS(windowIndex+1);
        windowSize = windowRight-windowLeft;

        % get data
        currentWindowIndices = ...
            find((myLengthSumNewborns{datasetIdx}>windowLeft) == (myLengthSumNewborns{datasetIdx}<windowRight));
        currentTotalLength = myLengthSumNewborns{datasetIdx}(currentWindowIndices);
        currentRatios = Ratios{datasetIdx}(currentWindowIndices);

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
        [counts,centers]=hist(currentRatios,mycenters);        
        modcounts = (windowSize*counts/sum(counts)+windowLeft);
        plot(modcounts,centers,'-','Color',myColors(windowIndex,:),'LineWidth',LINEWIDTH)
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

%% Plot lifetime against birth length

h=figure(5); clf; hold on
myColors = linspecer(numel(datasetsPaths));
myPlotMarkers = 'os^vd<>';

title([WHATDATA ' condition']);

%
meanInterDivisionTimes=[];
for datasetIdx = 1:numel(datasetsPaths)

    plot(allLengths{datasetIdx},lifeTimes{datasetIdx},'.',...
        'Marker',myPlotMarkers(datasetIdx),...
        'LineWidth',2,...
        'Color',myColors(datasetIdx,:));%,'MarkerFaceColor',myColors(datasetIdx,:));
    
    currentInterDivTimes = lifeTimes{datasetIdx};
    meanInterDivisionTimes(datasetIdx)=mean(currentInterDivTimes(~isnan(currentInterDivTimes)));            
        % note meanInterDivisionTimes is subject to how long cells where in
        % filamented state in this dataset.
    
end
%title(num2str(meanInterDivisionTimes))

% Plot the means per window (which was set in earlier plot)
for datasetIdx = 1:numel(datasetsPaths)
    
    %allSetsLifeTimes = [lifeTimes{:}];
    %allSetsallLengths = [allLengths{:}];
    
    meansLengths=[]; meansLifeTimes=[]; stdLifeTimes = [];
    for windowIndex = 1:(numel(WINDOWBORDERS)-1)

        windowLeft = WINDOWBORDERS(windowIndex);
        windowRight  = WINDOWBORDERS(windowIndex+1);
        windowSize = windowRight-windowLeft;

        % get data
        currentWindowIndices = ...
            find((allLengths{datasetIdx}>windowLeft) == (allLengths{datasetIdx}<windowRight));
        currentLengths = allLengths{datasetIdx}(currentWindowIndices);
        currentLifeTimes = lifeTimes{datasetIdx}(currentWindowIndices);

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
        'Marker',myPlotMarkers(datasetIdx),...
        'MarkerFaceColor','k','MarkerSize',10)
end
    
%xlim([LEFTX,RIGHTX]);
xlabel(['Birth size by ' LENGTHFIELD ' '],'Interpreter','None');
ylabel('Lifetime [min]');
MW_makeplotlookbetter(15)

if ~NOSAVEPLEASE
    saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_histograms.fig']);
    saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_histograms.tif']);
end

%% Plot with time-coded ratios
MINUTELIFETIMETRESHOLD=40;

h1=figure(6); clf; hold on
h2=figure(7); clf; hold on

% Get a colormap (is 64 long by default)
timeColormap = colormap(jet); % winter
rowsInColormap = size(timeColormap,1);

% 
longestTime = max([lifeTimes{:}]);
%longestTime = max([birthTimes{:}]);
fromLifeTimeToColorCodeFactor = longestTime/(rowsInColormap-1);

for datasetIdx = 1:numel(datasetsPaths)

    for i = 1:numel(myLengthSumNewborns{datasetIdx})

        %currentLifeTime = lifeTimes{datasetIdx}(myNewBornSchnitzNrs{datasetIdx}(i)); % daughter lifetime
        currentTime = lifeTimes{datasetIdx}(i); % parent lifetime
        %currentTime = birthTimes{datasetIdx}(i); % time at which born
        
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
        figure(6);
        plot(myLengthSumNewborns{datasetIdx}(i),Ratios{datasetIdx}(i),'o','MarkerSize',10,...
                'MarkerEdgeColor',colorForThisDataPoint,...
                'MarkerFaceColor',colorForThisDataPoint)
        
        % plot selection of datapoints to 2nd figure
        figure(7);
        if currentTime>MINUTELIFETIMETRESHOLD % || isnan(currentLifeTime)
            plot(myLengthSumNewborns{datasetIdx}(i),Ratios{datasetIdx}(i),'o','MarkerSize',10,...
                'MarkerEdgeColor',colorForThisDataPoint,...
                'MarkerFaceColor',colorForThisDataPoint)
        end
    end
    
end

for myfignum = 6:7
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

if ~NOSAVEPLEASE
    saveas(h1, [PLOTSAVEDIR WHATDATA '_ratiosByInterdiv.fig']);
    saveas(h1, [PLOTSAVEDIR WHATDATA '_ratiosByInterdiv.tif']);

    saveas(h2, [PLOTSAVEDIR WHATDATA '_ratiosByInterdivSelected.fig']);
    saveas(h2, [PLOTSAVEDIR WHATDATA '_ratiosByInterdivSelected.tif']);
end

%% Histogram of lifetimes
h=figure(8); clf; hold on

hist([lifeTimes{:}],30)
xlabel('Interdivision lifetime [min]');
ylabel('Count')

MW_makeplotlookbetter(20);

if ~NOSAVEPLEASE
    saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_generalhistogram.fig']);
    saveas(h, [PLOTSAVEDIR WHATDATA '_lifetime_generalhistogram.tif']);
end

%% histogram of cell sizes

figure(9); clf; hold on;
allLengths = [schnitzcells(:).(LENGTHFIELD)];
histogram(allLengths,100);

xlabel('Length (um)');
ylabel('Count');
title('Cell lenghts observed during movie');
MW_makeplotlookbetter(15);

%%
%{
figure(FIGURENUMBERS(end)+3); clf; hold on;

plot(histArea(2,:),histArea(1,:),'-r');
plot(histSkel(2,:),histSkel(1,:),'-b');
%}

%% 
disp('Done.');




