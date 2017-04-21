
% This script analyzes the Fts-labeled strain datasets.


%% param settings
WINDOWBORDERS   = [3:6:36];
BARCOLORS       = [240,155,34; 45 177 65; 37 156 190; 131 84 162; 241 88 58]./255;


%% Temperature datasets

if any(strcmp(RUNSECTIONSFILEFTS,{'all','loadData'}))

    % skeletonDataPaths{x} should match with schnitzDataPaths{x}.

    % stored skeleton analysis
    skeletonDataPaths={...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-07_asc777_temperatureRecovery\pos2crop\data\pos2crop-skeletonData.mat',...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-07_asc777_temperatureRecovery\pos3crop\data\pos3crop-skeletonData.mat'};
        % per frame the skeleton

    % stored fluorescence data
    fluorDataPaths={...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-07_asc777_temperatureRecovery\pos2crop\analysis\straightenedCells\2016-04-07pos2crop_straightFluorData.mat',...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-07_asc777_temperatureRecovery\pos3crop\analysis\straightenedCells\2016-04-07pos3crop_straightFluorData.mat'};
        %(note that kymograph was chosen from G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos2crop), first 20 frames
        %kymographs are made with script20160610_PetraFtsLocations    

    schnitzDataPaths={
        'G:\EXPERIMENTAL_DATA_2016\2016-04-07_asc777_temperatureRecovery\pos2crop\data\pos2crop-Schnitz.mat',...
        'G:\EXPERIMENTAL_DATA_2016\2016-04-07_asc777_temperatureRecovery\pos3crop\data\pos3crop-Schnitz.mat'};
end    
    
%% convert schnitzcells data to ages per frame per cell
if any(strcmp(RUNSECTIONSFILEFTS,{'all','getAges'}))

    allAgesPerFrame = {};
    for datasetIndex = 1:numel(schnitzDataPaths)

        load(schnitzDataPaths{datasetIndex});
        allAgesPerFrame{datasetIndex} = {};
        for schnitzIdx = 1:numel(schnitzcells)
            for frameIdx = 1:numel(schnitzcells(schnitzIdx).frame_nrs)

                % get frame nr and cell nr.
                currentFrameNr = schnitzcells(schnitzIdx).frame_nrs(frameIdx); 
                currentCellNo  = schnitzcells(schnitzIdx).cellno(frameIdx);

                % store age of cell accordingly
                allAgesPerFrame{datasetIndex}{currentFrameNr}(currentCellNo) = schnitzcells(schnitzIdx).time(frameIdx) - schnitzcells(schnitzIdx).birthTime;

            end
        end
    
    end
        
end

%%

if any(strcmp(RUNSECTIONSFILEFTS,{'all','analyzeData'}))

    PEAKTRESHOLD=400;
    figure(1); clf;
    figure(101); clf;
    plotcolors = 'rrr';
    scatterX={}; scatterY={}; 
    lifeTimesScatter = {};

    PLOTHELPINGLINES=1;

    %% Collect all peak locations

    for datasetIndex = 1:numel(skeletonDataPaths)

        %%
        %load ([p.tracksDir p.movieName '-skeletonData.mat']);
        clear allpeakmeanY allpeakY;
        load (skeletonDataPaths{datasetIndex});

        %%
        %load(saveLocationMatFile([p.analysisDir 'straightenedPlots\' p.movieDate p.movieName '_straightFluorData.mat']);
        load (fluorDataPaths{datasetIndex});


        %%

        figure(1); hold on;
        %clf;
        %hold on;

        pileofXmicron=[];
        pileofXPx = [];
        pileofmeanY=[];
        selectionVector=[];
        pileofLengthsOfBacteriaInMicrons=[];
        correspondingLengthsToLocations=[];
        correspondingLengthsToLocationsPx=[];
        correspondingLifeTimesToLocations=[];
        for framenr=1:numel(allpeakXMicrons)

            if exist('allpeakY','var') % oops, have to resolve this later.. double param naming.. TODO
                meanYFrame  = allpeakY{framenr};
            else
                meanYFrame  = allpeakmeanY{framenr};
            end

            xpeaksFrame  = allpeakXMicrons{framenr};
            lengthsFrame = allLengthsOfBacteriaInMicrons{framenr};

            xpeaksFramePx  = allpeakXPixels{framenr};
            lengthsFramePx = allLengthsOfBacteriaInPixels{framenr};

            %disp(['adding ' num2str(numel(xpeaksFrame)) ' cells for frame ' num2str(framenr) '..']);
            pileofLengthsOfBacteriaInMicrons =...
                [pileofLengthsOfBacteriaInMicrons lengthsFrame];

            for cellno = 1:numel(xpeaksFrame)

                %plot(xpeaksFrame)
                if ~(isempty(xpeaksFrame) || isempty(meanYFrame))

                    % raw data

                    sizexpeaksParam = size(xpeaksFrame{cellno});
                    if sizexpeaksParam(1) > sizexpeaksParam(2) % oops, have to resolve this later.. double param standards.. TODO
                        pileofXmicron = [pileofXmicron xpeaksFrame{cellno}'];
                        pileofXPx     = [pileofXPx xpeaksFramePx{cellno}'];
                    else
                        pileofXmicron = [pileofXmicron xpeaksFrame{cellno}];
                        pileofXPx     = [pileofXPx xpeaksFramePx{cellno}];
                    end

                    fluorPeaksThisCell = meanYFrame{cellno};
                    pileofmeanY = [pileofmeanY fluorPeaksThisCell];

                    % create duplicate lengths for datapoints belonging to same
                    % bacteria
                    duplicatedLengths = ones(1,numel(xpeaksFrame{cellno}))*lengthsFrame(cellno);
                    correspondingLengthsToLocations = [correspondingLengthsToLocations duplicatedLengths];
                    % same w. pixels
                    duplicatedLengthsPx = ones(1,numel(xpeaksFrame{cellno}))*lengthsFramePx(cellno);
                    correspondingLengthsToLocationsPx = [correspondingLengthsToLocationsPx duplicatedLengthsPx];
                    % same w. lifetime cell
                    duplicatedLifetimes = ones(1,numel(xpeaksFrame{cellno}))*allAgesPerFrame{datasetIndex}{framenr}(cellno);
                    correspondingLifeTimesToLocations = [correspondingLifeTimesToLocations duplicatedLifetimes];
                    
                    % select peaks that are >200
                    selectionVector = [selectionVector fluorPeaksThisCell>PEAKTRESHOLD];

                end
                %selectedataXmicron = ...

            end

        end

        %plot(correspondingLengthsToLocations,pileofXmicron,'.')
        %plot(correspondingLengthsToLocations,pileofXmicron./correspondingLengthsToLocations,'.')
        indicesToSelect = find(selectionVector);
        scatterX{datasetIndex}=correspondingLengthsToLocations(indicesToSelect);
        scatterY{datasetIndex}=pileofXmicron(indicesToSelect)./correspondingLengthsToLocations(indicesToSelect);
        lifeTimesScatter{datasetIndex} = correspondingLifeTimesToLocations(indicesToSelect);
        plot(scatterX{datasetIndex},scatterY{datasetIndex},'x','Color',plotcolors(datasetIndex))

        % and in Pixels
        figure(101); hold on;
        scatterXPx{datasetIndex}=correspondingLengthsToLocationsPx(indicesToSelect);
        scatterYPx{datasetIndex}=pileofXPx(indicesToSelect)./correspondingLengthsToLocationsPx(indicesToSelect);
        plot(scatterXPx{datasetIndex},scatterYPx{datasetIndex},'.','Color',plotcolors(datasetIndex))


        %plot(allpeak)

    end
    
end

%% Sanity check
if any(strcmp(RUNSECTIONSFILEFTS,{'all','LengthVsDivisionLocation'}))

    for datasetIndex = 1:numel(skeletonDataPaths)
    
        MICRONSPERPIXEL=0.0431;

        figure(1); hold on;
        plot(scatterXPx{datasetIndex}*MICRONSPERPIXEL,scatterYPx{datasetIndex},'.b')

    end
    clear datasetIndex
    %% figure 2, make a scatter plot (ratio S vs length L)

    MARKERCOLOR = [.5 .5 .5];
    ALPHA = .10;
    NRCONTOURLINES=0;
    MAKEUSEOFSYMMETRY=1;
    COLOR='b';

    h=figure(2); clf; hold on;
    set(h,'Position',[200,200,600+200,400+200]);

    if ~MAKEUSEOFSYMMETRY
        data = [[scatterX{:}]; [scatterY{:}]]';
    else
        data = [[[scatterX{:}] [scatterX{:}]]; [[scatterY{:}], 1-[scatterY{:}]]]'; % symmetry data
    end
       
    if ~exist('PAINTWITHTIME','var')
        scatter([scatterX{:}],[scatterY{:}],12^2,'filled','MarkerFaceColor',MARKERCOLOR,'MarkerFaceAlpha',ALPHA);
        %scatter([scatterX{:}],[scatterY{:}],75,'filled','MarkerFaceColor',[.7 .7 .7]);%,'MarkerFaceAlpha',3/8);
        if MAKEUSEOFSYMMETRY
            scatter([scatterX{:}],1-[scatterY{:}],12^2,'filled','MarkerFaceColor',MARKERCOLOR,'MarkerFaceAlpha',ALPHA); % plot symmetric
            %scatter([scatterX{:}],1-[scatterY{:}],75,'filled','MarkerFaceColor',[.7 .7 .7]); % plot symmetric
        end
    else
        % calcalate ages matching colors
        % Note that a load of data comes from cells that are growing but not
        % dividing, so the age of these cells is meaningless; rearrangements
        % are due to growth..
        allTimeDataScatter = [lifeTimesScatter{:}];
        maxAllTimeDataScatter=max(allTimeDataScatter);
        allTimeDataScatterNormalized = allTimeDataScatter./maxAllTimeDataScatter;
        myTimecolormap = [makeColorMap([1 0 0],[0 1 0],[0 0 1]); repmat([0 0 1],200,1)];
        allTimesColors = myTimecolormap(int64(allTimeDataScatterNormalized*299)+1,:);


        scatter([scatterX{:}],[scatterY{:}],12^2,allTimesColors);%'filled','MarkerFaceColor',[.7 .7 .7],'MarkerFaceAlpha',0.15);
        colormap(myTimecolormap);  

        hcb=colorbar;

        inputSettings.rangeIn = [0,1];
        inputSettings.desiredSpacing = 50;
        inputSettings.rangeOut = [0,maxAllTimeDataScatter];
        [tickLocationsOldMetric, correspondingLabels] = labelremapping(inputSettings);

        set(hcb,'YTick',tickLocationsOldMetric,'YTickLabel',correspondingLabels)    
        
        %{
        % Plot cells of certain age (young) -- comment out above scatter to use
        youngCellIndices = allTimeDataScatter<90;
        xData=[scatterX{:}];
        yData=[scatterY{:}];
        scatter(xData(youngCellIndices),yData(youngCellIndices),12^2,allTimesColors(youngCellIndices));%'filled','MarkerFaceColor',[.7 .7 .7],'MarkerFaceAlpha',0.15);
        scatter(xData(youngCellIndices),1-yData(youngCellIndices),12^2,allTimesColors(youngCellIndices));%'filled','MarkerFaceColor',[.7 .7 .7],'MarkerFaceAlpha',0.15);
        %}
    end
    %}

    %{
    plot([scatterX{:}],[scatterY{:}],['.' COLOR],'MarkerSize',7);
    if MAKEUSEOFSYMMETRY
        plot([scatterX{:}],1-[scatterY{:}],['.' COLOR],'MarkerSize',7); % plot symmetric
    end
    %}

    [bandwidth,density,X,Y] = kde2d(data);      
    %[C, l1] = contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);

    % density corrected for data loss at higher lengths..
    %[C, l1] = contour(X,Y,density.*X,NRCONTOURLINES,'-k','LineWidth',2);

    if ~exist('DONTPLOTCONTOURLINES','var')
        [C, l1] = contour(X,Y,density.*X,NRCONTOURLINES,'-k','LineWidth',2);
    end

    %title('kde, density directed (pdf(x)*x)')
    xlabel(['Length of cell [' '\mu' 'm]']);
    ylabel(['Relative FtsA peak locations']);
    %figure(); hist([scatterX{:}])

    maxX = max([scatterX{:}]);
    xlim([0, maxX]);
    ylim([0,1]);

    MW_makeplotlookbetter(24);

    % plot helping lines at 1/2n
    %{
    if PLOTHELPINGLINES
        N=5;
        for windowIndex=1:N
            for j = 1:(windowIndex*2-1)
                plot([0, maxX], [(j)/(2*windowIndex) (j)/(2*windowIndex)],'-','Color',[.5 .5 .5],'LineWidth',N-windowIndex+1)
            end
        end
    end
    %}

    % plot helping lines at 1/2n
    if PLOTHELPINGLINES
        N=5;
        for windowIndex=1:(numel(WINDOWBORDERS)-1)
            for j = 1:2:(windowIndex*2-1)
                plot([WINDOWBORDERS(windowIndex), WINDOWBORDERS(windowIndex+1)], [(j)/(2*windowIndex) (j)/(2*windowIndex)],'-','Color',BARCOLORS(windowIndex,:),'LineWidth',8)
            end
        end
    end

    if exist('CUSTOMXLIM','var');
       xlim(CUSTOMXLIM);
    end

    SIZE=[7.5 6.5];
    OFFSET = [3 3];
    h.Units = 'centimeters';
    h.Position = [OFFSET SIZE]*2;
    h.PaperUnits = 'centimeters';
    h.PaperPosition = [0 0 SIZE]*2;


    %OUTPUTDIRfts='D:\Local_Data\Dropbox\Dropbox\Filamentation recovery\MW\figures_new\matlab_export\FtsaSPlot\';
    if exist('OUTPUTDIRfts','var') 
        saveas(h,[OUTPUTDIRfts 'rutger_ftsa_peaks_new.tif'])
        saveas(h,[OUTPUTDIRfts 'rutger_ftsa_peaks_new.svg'])
        saveas(h,[OUTPUTDIRfts 'rutger_ftsa_peaks_new.fig'])
    else
        disp('Set OUTPUTDIRfts to save plot');
    end
    % saveas(2, 'D:\Local_Data\Dropbox\Dropbox\Filamentation recovery\MW\figures_new\matlab_export\rutger_ftsa_peaks_new.svg')

    hFig3c = h;
    
end
    
%% Create overlay plot if division ratios are avaible

if any(strcmp(RUNSECTIONSFILEFTS,{'all','divisionLocations'}))

    NEWMAXX = 50;

    figure(3); clf; hold on;
    [C, l1] = contour(X,Y,density.*X,NRCONTOURLINES,'-k','LineWidth',2);

    % data that should be obtained from
    % script20160429_filamentRecoveryDivisionRatioss
    if exist('Ratios','var')
        for dataSetIndex = 1:numel(datasetsPaths)
            plot(myLengthSumNewborns{dataSetIndex},Ratios{dataSetIndex},'o', 'Color', PLOTCOLORS(dataSetIndex,:),'LineWidth',2);
        end
    end

    xlim([0, NEWMAXX]);
    ylim([0,1]);

    xlabel(['Length of cell [' '\mu' 'm]']);
    ylabel('ftsA peaks location / division location');

    MW_makeplotlookbetter(15)
    
end

%% length histograms, also per regime
if any(strcmp(RUNSECTIONSFILEFTS,{'all','lengthHistogramAndPerGroup'}))

    figure(4);clf;

    [n,c] = hist(pileofLengthsOfBacteriaInMicrons,50);
    plot(c,n./(sum(n)*(c(2)-c(1))),'o-k','MarkerSize',15,'LineWidth',3);
    MW_makeplotlookbetter(20);
    xlabel('Cell length [um]');
    ylabel('Probability');

    %% ring statistics (locations overrepresented?)

    NUMBINS=20;

    % 14-23 um regime
    % scatterX{2}, scatterY{2} give best current available data

    dataSetIndex=2;

    regionsDouble = {[0,0.5],[0,0.5],[0,2/6,.5],[0,2/8,.5],[0,2/10,4/10,5/10]};

    theBinEdges = [0:.5/NUMBINS:0.5];
    thebinCenters=[theBinEdges(2:end)+theBinEdges(1:end-1)]./2;

    %WINDOWBORDERS = [14,23,31];
    %regions = {[0,2/6,4/6,1],[0,2/8,4/8,6/8,8/8]};
    %regionsDouble = {[0,2/6,.5],[0,2/8,.5]};

    eventsPerRatio = {}; pdf={}; normalizedPdf=  {}; ftsData=struct;
    for dataSetIndex = 1:numel(skeletonDataPaths)
        for windowIndex=1:(numel(WINDOWBORDERS)-1)

            relevantDataIndices = (scatterX{dataSetIndex}>WINDOWBORDERS(windowIndex) & scatterX{dataSetIndex}<WINDOWBORDERS(windowIndex+1));
            selectedXdata = scatterX{dataSetIndex}(relevantDataIndices);
            selectedYdata = scatterY{dataSetIndex}(relevantDataIndices);

            % make use of symmetery
            doubleDataX = [selectedXdata 1-selectedXdata];
            doubleDataY = [selectedYdata 1-selectedYdata];

            %plot(selectedXdata,selectedYdata,'.');
            % how many divisions occured at one of the ratios in this regime?
            % E.g . how many events occured inbetween 0 and 2/6 for regime 3?
            eventsPerRatio = histcounts(doubleDataY,regionsDouble{windowIndex}); % note histogram() and histcounts() take edges as input, hist() does not

            [count, theBinEdges] = histcounts(doubleDataY,theBinEdges); % note histcounts() takes edges as input, hist() does not
            %count=h.Values;

            dt=theBinEdges(2)-theBinEdges(1);
            normalizedpdf = count./sum(count)*dt;

            %plot(thebinCenters,count,'-','LineWidth',4);

            disp(['Stats for regime ' num2str(windowIndex) ': ' num2str(WINDOWBORDERS(windowIndex)) 'um -' num2str(WINDOWBORDERS(windowIndex+1)) 'um' ]);
            disp(num2str(eventsPerRatio));

            ftsData.count{dataSetIndex}{windowIndex}           = count; % 
            ftsData.eventsPerRatio{dataSetIndex}{windowIndex} = eventsPerRatio; % 
            ftsData.normalizedPdf{dataSetIndex}{windowIndex}   = normalizedPdf;

        end
    end

    ftsData.centers         = thebinCenters;

    % Store data and also create summary (mean, sum, normalized) params and store those
    for windowIndex = 1:(numel(WINDOWBORDERS)-1)
        % determine average function
        countsFromMultipleDatasets = arrayfun(@(x) ftsData.count{x}{windowIndex}, 1:numel(ftsData.count), 'UniformOutput', false);
        ftsData.meanCounts{windowIndex}=mean(cell2mat(countsFromMultipleDatasets'));
        ftsData.sumCounts{windowIndex}=sum(cell2mat(countsFromMultipleDatasets'));

        % normalize the pdf
        dt=ftsData.centers(2)-ftsData.centers(1);
        ftsData.normalizedPdf{windowIndex} = ftsData.sumCounts{windowIndex}./sum(ftsData.sumCounts{windowIndex})*dt;

        % sum the different event counts for the ratios
        ratiocountsFromMultipleDatasets = arrayfun(@(x) ftsData.eventsPerRatio{x}{windowIndex}, 1:numel(ftsData.eventsPerRatio), 'UniformOutput', false);
        ftsData.sumRatioCounts{windowIndex}=sum(cell2mat(ratiocountsFromMultipleDatasets'));
    end

    % plot it
    figure(5); clf; hold on;

    if ~exist('toPlotWindows')
        toPlotWindows = 1:(numel(WINDOWBORDERS)-1);
    end
    for windowIndex=toPlotWindows
        plot(ftsData.centers,ftsData.normalizedPdf{windowIndex},'-','LineWidth',4);
    end

    xlabel('Relative ring location, S_{ring}');
    ylabel(['Probability']);

    xlim([0,0.5])

    %ylim([0,max([pdf{:}])*1.2]);

    MW_makeplotlookbetter(15);
    
end

%% 
% Execute
% script20160429_filamentRecoveryDivisionRatioss
% first

if any(strcmp(RUNSECTIONSFILEFTS,{'all','compareDivisionsWRings'}))

    REGIME=3;

    colors=linspecer(2);
    figure(6); clf; hold on;

    theYlim = [0, max([histData.normalizedPdf{REGIME}*2 ftsData.normalizedPdf{REGIME}])*1.2];

    % plot ratios
    for fractions = 1:2:(2*REGIME)
        ratio = fractions/(REGIME*2);
        plot([ratio,ratio], theYlim,'-','Color',[.7 .7 .7],'LineWidth',3)
    end

    % observed rings
    l2=plot(ftsData.centers,ftsData.normalizedPdf{REGIME},'-','Color',colors(1,:),'LineWidth',3);
    %bar(ftsData.centers, ftsData.normalizedPdf{REGIME}),'FaceColor',[.7 .7 .7]);
    % observed divisions 
    l1=plot(histData.centers, histData.normalizedPdf{REGIME}*2,'-','Color','k','LineWidth',3); % 2* since renormalization to new xlimits
    %bar(histData.centers, histData.normalizedPdf{REGIME}*2,'FaceColor','k');

    % cosmetics
    legend([l1,l2],{'Observed divisions','Observed rings'});
    MW_makeplotlookbetter(15);
    xlim([0,.5]);
    ylim(theYlim)

    xlabel('Relative division location, S')
    %xlabel('L_{daughter}/L_{mother}');
    ylabel('Probability');

    % prints stats
    disp(['Div/FtsA stats for regime ' num2str(REGIME) ': ' num2str(WINDOWBORDERS(REGIME)) 'um -' num2str(WINDOWBORDERS(REGIME+1)) 'um' ]);
    disp(num2str(histData.sumRatioCounts{REGIME}));
    disp(num2str(ftsData.sumRatioCounts{REGIME}));

end
%% more sanity checks

%load (['G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos' num2str(datasetIndex) 'crop\data\pos' num2str(datasetIndex) 'crop-Schnitz.mat']);







