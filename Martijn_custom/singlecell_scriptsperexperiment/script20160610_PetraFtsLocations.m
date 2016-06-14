



%%
PEAKTRESHOLD=400;
figure(1); clf;
figure(101); clf;
plotcolors = 'rrr';
scatterX={}; scatterY={}; 


%%
for datasetIndex = 1

    %%
    %load ([p.tracksDir p.movieName '-skeletonData.mat']);
    load (['G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-19_part2\pos' num2str(datasetIndex) 'crop\data\pos' num2str(datasetIndex) 'crop-skeletonData.mat']);

    %%
    %load(saveLocationMatFile([p.analysisDir 'straightenedPlots\' p.movieDate p.movieName '_straightFluorData.mat']);    
    load(['G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-19_part2\pos' num2str(datasetIndex) 'crop\analysis\straightenedCells\2016-05-19pos' num2str(datasetIndex) 'crop_straightFluorData.mat']);


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
    for framenr=1:numel(allpeakXMicrons)

        meanYFrame  = allpeakmeanY{framenr};
        xpeaksFrame  = allpeakXMicrons{framenr};
        lengthsFrame = allLengthsOfBacteriaInMicrons{framenr};

        xpeaksFramePx  = allpeakXPixels{framenr};
        lengthsFramePx = allLengthsOfBacteriaInPixels{framenr};
        
        %disp(['adding ' num2str(numel(xpeaksFrame)) ' cells for frame ' num2str(framenr) '..']);
        pileofLengthsOfBacteriaInMicrons =...
            [pileofLengthsOfBacteriaInMicrons lengthsFrame];

        for cellno = 1:numel(xpeaksFrame)

            %plot(xpeaksFrame)
            if ~isempty(xpeaksFrame)

                % raw data
                pileofXmicron = [pileofXmicron xpeaksFrame{cellno}'];
                pileofXPx     = [pileofXPx xpeaksFramePx{cellno}'];
                fluorPeaksThisCell = meanYFrame{cellno};
                pileofmeanY = [pileofmeanY fluorPeaksThisCell];

                % create duplicate lengths for datapoints belonging to same
                % bacteria
                duplicatedLengths = ones(1,numel(xpeaksFrame{cellno}))*lengthsFrame(cellno);
                correspondingLengthsToLocations = [correspondingLengthsToLocations duplicatedLengths];
                % same w. pixels
                duplicatedLengthsPx = ones(1,numel(xpeaksFrame{cellno}))*lengthsFramePx(cellno);
                correspondingLengthsToLocationsPx = [correspondingLengthsToLocationsPx duplicatedLengthsPx];
                
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
    plot(scatterX{datasetIndex},scatterY{datasetIndex},'x','Color',plotcolors(datasetIndex))

    % and in Pixels
    figure(101); hold on;
    scatterXPx{datasetIndex}=correspondingLengthsToLocationsPx(indicesToSelect);
    scatterYPx{datasetIndex}=pileofXPx(indicesToSelect)./correspondingLengthsToLocationsPx(indicesToSelect);
    plot(scatterXPx{datasetIndex},scatterYPx{datasetIndex},'.','Color',plotcolors(datasetIndex))
    
    
    %plot(allpeak)

end

%% Sanity check
MICRONSPERPIXEL=0.0431;

figure(1); hold on;
plot(scatterXPx{datasetIndex}*MICRONSPERPIXEL,scatterYPx{datasetIndex},'.b')

%% 
NRCONTOURLINES=3;
MAKEUSEOFSYMMETRY=1;
COLOR='b';

figure(2); clf; hold on;

if ~MAKEUSEOFSYMMETRY
    data = [[scatterX{:}]; [scatterY{:}]]';
else
    data = [[[scatterX{:}] [scatterX{:}]]; [[scatterY{:}], 1-[scatterY{:}]]]'; % symmetry data
end

plot([scatterX{:}],[scatterY{:}],['.' COLOR],'MarkerSize',7);
if MAKEUSEOFSYMMETRY
    plot([scatterX{:}],1-[scatterY{:}],['.' COLOR],'MarkerSize',7); % plot symmetric
end
[bandwidth,density,X,Y] = kde2d(data);      
%[C, l1] = contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);

% density corrected for data loss at higher lengths..
%[C, l1] = contour(X,Y,density.*X,NRCONTOURLINES,'-k','LineWidth',2);
[C, l1] = contour(X,Y,density.*X,NRCONTOURLINES,'-k','LineWidth',2);

%title('kde, density directed (pdf(x)*x)')
xlabel(['Length of cell [' '\mu' 'm]']);
ylabel('ftsA peaks location');
%figure(); hist([scatterX{:}])

maxX = max([scatterX{:}]);
xlim([0, maxX]);
ylim([0,1]);

MW_makeplotlookbetter(20);

% plot helping lines at 1/2n
N=5;
for i=1:N
    for j = 1:(i*2-1)
        plot([0, maxX], [(j)/(2*i) (j)/(2*i)],'-','Color',[.5 .5 .5],'LineWidth',N-i+1)
    end
end

%% Create overlay plot if division ratios are avaible
NEWMAXX = 50;

figure(3); clf; hold on;
[C, l1] = contour(X,Y,density.*X,NRCONTOURLINES,'-k','LineWidth',2);

% data that should be obtained from
% script20160429_filamentRecoveryDivisionRatioss
if exist('Ratios','var')
    for datasetIdx = 1:numel(datasetsPaths)
        plot(myLengthSumNewborns{datasetIdx},Ratios{datasetIdx},'o', 'Color', PLOTCOLORS(datasetIdx,:),'LineWidth',2);
    end
end

xlim([0, NEWMAXX]);
ylim([0,1]);

xlabel(['Length of cell [' '\mu' 'm]']);
ylabel('ftsA peaks location / division location');

MW_makeplotlookbetter(15)

%% 

figure(4);clf;

[n,c] = hist(pileofLengthsOfBacteriaInMicrons,50);
plot(c,n./(sum(n)*(c(2)-c(1))),'o-k','MarkerSize',15,'LineWidth',3);
MW_makeplotlookbetter(20);
xlabel('Cell length [um]');
ylabel('Probability');

%% Create kymograph -------------------------------------------------------

%% Create a lookup table
%% build a lookuptable that corresponds with lineage structure
% lookuptable{n} gives, for nth simulation time, in lookuptable{n}(:,1) the
% schnitzes that are alive during that time, and lookuptable{n}(:,2), the
% corresponding point in their cellphase.
% lookuptable{n} = [schnitznr, framenr]

alltimes = unique([schnitzcells.time]);

lookuptable={}; lookuptable{1}=[]; hits=0;
% for first frame
for i = 1:numel(schnitzcells)
    if any(schnitzcells(i).time==alltimes(1))
        hits=hits+1;
        lookuptable{1}(hits,1) = i;
        lookuptable{1}(hits,2) = 1;
    end
end

tic;
% continue with all consequtive frames
for frameNr = 2:numel(alltimes)
    
    lookuptable{frameNr}=[]; 
    rprime = 0;
    for r = 1:numel(lookuptable{frameNr-1}(:,1))
        
        % previous schnitz
        previousSchnitz = lookuptable{frameNr-1}(r,1);
        
        % previous indices of time
        previousTimeIdx = find(schnitzcells(previousSchnitz).time==alltimes(frameNr-1));
        
        % if previous time is last time
        if numel((schnitzcells(previousSchnitz).time))==previousTimeIdx
            % cell has divided
            nextSchnitz1 = schnitzcells(previousSchnitz).D;
            nextSchnitz2 = schnitzcells(previousSchnitz).E;
            % put in table
            rprime=rprime+1;
            lookuptable{frameNr}(rprime,1) = nextSchnitz1;
            lookuptable{frameNr}(rprime,2) = 1;
            rprime=rprime+1;
            lookuptable{frameNr}(rprime,1) = nextSchnitz2;
            lookuptable{frameNr}(rprime,2) = 1;
        else
            % cell has not divided
            nextSchnitz = previousSchnitz;
            rprime=rprime+1;
            lookuptable{frameNr}(rprime,1) = nextSchnitz;
            lookuptable{frameNr}(rprime,2) = find((schnitzcells(nextSchnitz).time==alltimes(frameNr)));
        end
            
        if toc>10
            error('timeout error');
        end
        
    end
    
end

%% Now plot the kymograph
LENGTHFIELD = 'pixLength_skeleton';
RANGE = 1:100;
MAXPIXELHEIGHTKYMOGRAPH=500;

% determine how many timepoints need to be plotted
cellPhaseActiveIdxs=[];
for frameNr = RANGE
    if ~isempty(allmeanY{frameNr})
        cellPhaseActiveIdxs(end+1) = frameNr;
    end
end
sizeX = numel(cellPhaseActiveIdxs); % frames 

% determine total length to be plotted
theLastFrameNr = max(cellPhaseActiveIdxs);
lastSchnitzNrs = lookuptable{theLastFrameNr}(:,1)';
summedNumel = 0; allYnumels = [];
for schnitzIdx=1:numel(lastSchnitzNrs)
    % currentLength = schnitzcells(lastSchnitzNrs).(LENGTHFIELD);
        
    currentCellPhaseIdx = lookuptable{theLastFrameNr}(schnitzIdx,2);
    
    cellno=schnitzcells(lastSchnitzNrs(schnitzIdx)).cellno(currentCellPhaseIdx);
    currentNumelYvalues = numel(allmeanY{theLastFrameNr}{cellno});
    
    allYnumels(end+1)=currentNumelYvalues;
    summedNumel = summedNumel + currentNumelYvalues + 1; % +1 to allow for separator pixel
end

sizeY = min(MAXPIXELHEIGHTKYMOGRAPH,uint16(summedNumel)); % cell length
outputMatrix = zeros(sizeX,sizeY,3);%% 


% dummy markers for legend (appear at t<0)
%l1=plot(-1,0,'sk-','MarkerFaceColor','k','MarkerSize',5);
%l2=plot(-1,0,'or','MarkerSize',4,'MarkerFaceColor','r');
%l3=plot(-1,0,'ob','MarkerSize',4);

%imshow(outputMatrix',[]);

%% normalize Y vector (fluor vector)

allYdata1 = [allmeanY{:}];
allYdata2 = [allYdata1{:}];
minYdata=min(allYdata2);
maxYdata=max(allYdata2);
rangeYdata = maxYdata-minYdata;

normalizedAllmeanY = cell(size(allmeanY));
for i=1:numel(allmeanY)
    for j = 1:numel(allmeanY{i})
        normalizedAllmeanY{i}{j} = (allmeanY{i}{j}-minYdata) ./ rangeYdata;
    end
end

workAllmeanY = normalizedAllmeanY;

%% create graphical representation plot & gather info for kymograph
figure(5); clf; hold on;

tic
for frameIdx = 1:numel(cellPhaseActiveIdxs)
    
    currentSchnitzes    = lookuptable{cellPhaseActiveIdxs(frameIdx)}(:,1)';
    currentCellPhaseIdx = lookuptable{cellPhaseActiveIdxs(frameIdx)}(:,2)';
    
    totalLength=0;
    for loopIdx = 1:numel(currentSchnitzes)
        totalLength = totalLength+schnitzcells(currentSchnitzes(loopIdx)).(LENGTHFIELD)(currentCellPhaseIdx(loopIdx));
    end
    
    previousLengthsSummed=0;
    previousRawPixelLengthsSummed=0;
    plottingDividedLocations = [-totalLength/2];
    plottingRingLocations = []; plottingNucleoidLocations= [];
    nrSeparatorPixels = 0;
    for loopIdx = 1:numel(currentSchnitzes)
    
        plottingDividedLocations = ...
            [plottingDividedLocations previousLengthsSummed+schnitzcells(currentSchnitzes(loopIdx)).(LENGTHFIELD)(currentCellPhaseIdx(loopIdx))-totalLength/2];
        
        cellno=schnitzcells(currentSchnitzes(loopIdx)).cellno(currentCellPhaseIdx(loopIdx));
        Yvector = workAllmeanY{cellPhaseActiveIdxs(frameIdx)}{cellno};
        currentRawPixelLength = numel(Yvector);
        for pxIdx = 1:currentRawPixelLength
            matrixIndex = pxIdx+previousRawPixelLengthsSummed+nrSeparatorPixels;
            % don't draw outside specified size
            if (matrixIndex)>MAXPIXELHEIGHTKYMOGRAPH
                break;
            end
            outputMatrix(frameIdx,matrixIndex,1) = Yvector(pxIdx);
            outputMatrix(frameIdx,matrixIndex,2) = Yvector(pxIdx);
            outputMatrix(frameIdx,matrixIndex,3) = Yvector(pxIdx);
        end
        if ~((matrixIndex)>MAXPIXELHEIGHTKYMOGRAPH)
            nrSeparatorPixels = nrSeparatorPixels+1;
            if currentCellPhaseIdx(loopIdx)==1
                outputMatrix(frameIdx,matrixIndex+1,1) = 1;
            else
                outputMatrix(frameIdx,matrixIndex+1,3) = 1;
            end
        end        
        
        %{
        for ringIdx = 1:numel(schnitzcells(currentSchnitzes(loopIdx)).absoluteRingSite{currentFrames(loopIdx)})
            plottingRingLocations = ...
                [plottingRingLocations previousLengthsSummed+schnitzcells(currentSchnitzes(loopIdx)).absoluteRingSite{currentFrames(loopIdx)}(ringIdx)-totalLength/2];
        end

        for nucleoidIdx = 1:numel(schnitzcells(currentSchnitzes(loopIdx)).absoluteNucleoidSite{currentFrames(loopIdx)})
            plottingNucleoidLocations = ...
                [plottingNucleoidLocations previousLengthsSummed+schnitzcells(currentSchnitzes(loopIdx)).absoluteNucleoidSite{currentFrames(loopIdx)}(nucleoidIdx)-totalLength/2];
        end
        %}
        
        %plot(schnitzcells(currentSchnitzes(loopIdx)).times(currentFrames(loopIdx)),...
        %     previousLengthsSummed+schnitzcells(currentSchnitzes(loopIdx)).(LENGTHFIELD)(currentFrames(loopIdx))-totalLength/2,...
        %    'ok-','MarkerFaceColor','k','MarkerSize',5);
        
        previousLengthsSummed=previousLengthsSummed+schnitzcells(currentSchnitzes(loopIdx)).(LENGTHFIELD)(currentCellPhaseIdx(loopIdx));
        previousRawPixelLengthsSummed=previousRawPixelLengthsSummed+currentRawPixelLength;
       
        if toc>10
            error('timeout error');
        end
        
    end
    
    plot(ones(1,numel(plottingDividedLocations))*alltimes(cellPhaseActiveIdxs(frameIdx)),...
            plottingDividedLocations,...
           'sk-','MarkerFaceColor','k','MarkerSize',5);
       %{
    plot(ones(1,numel(plottingRingLocations))*simulationtimes(t),...
            plottingRingLocations,...
           'or','MarkerSize',4,'MarkerFaceColor','r');
    if SHOWNUCLEOIDS
        plot(ones(1,numel(plottingNucleoidLocations))*simulationtimes(t),...
                plottingNucleoidLocations,...
               'ob','MarkerSize',4);
           %}
end

disp('done');

%% plot kymograph

figure(6); clf; hold on;
warning('Cell length might be distorted due to diagonal pixels..');
imshow(outputMatrix,[]);

imwrite(outputMatrix,'G:\kymograph.tif')

%% Create a little kymograph per schnitz

theoutputdir = [p.analysisDir '\straightenedCells\kymospercell\'];
if ~exist(theoutputdir,'dir')
    mkdir(theoutputdir);
end

% prepare ring timing plot
mycolor = linspecer(1);

multipleringTimePlotx={};
multipleringTimePloty={};

for schnitzIdx = 1:numel(schnitzcells)
    %%
   
    % Don't consider schnitzes that have unfinished life cycle
    if schnitzcells(schnitzIdx).E==0
        disp(['Skipping barren schnitzcell ' num2str(schnitzIdx)]);
        continue;
    end
    
    % prepare ring distribution plot
    h8=figure(8); clf; hold on;
    xlabel('Biased cell length [pixels]');
    ylabel('Fluor intensity (normalized by overall mean)');
    title('Grey towards black is cell cycle progression.');
    MW_makeplotlookbetter(15);
    
    % for ring timing plots
    ringTimePlotx=[];    
    ringTimePloty=[];
    
    % frame and cellnos applicable for this schnitz needed for calculations
    theframes  = schnitzcells(schnitzIdx).frame_nrs;
    thecellnos = schnitzcells(schnitzIdx).cellno;
    
    cellLifeLength = numel(theframes); % in frames, for administration
    
    allCellPixelLengths=[]; cellPhaseActiveIdxs = [];
    for i = 1:numel(theframes)
        if~isempty(allmeanY{theframes(i)})
            allCellPixelLengths(end+1) = numel(allmeanY{theframes(i)}{thecellnos(i)});
            cellPhaseActiveIdxs(end+1) = i;
        end
    end
    
    maxPixLength = max(allCellPixelLengths);
    nrActiveCellPhaseFrames = numel(cellPhaseActiveIdxs);
    
    littleOutputMatrix = zeros(nrActiveCellPhaseFrames,maxPixLength,3); %% 
    littleOutputMatrix(:,:,3) = ones(nrActiveCellPhaseFrames,maxPixLength);
        
    for idx = 1:nrActiveCellPhaseFrames        
        
        cellPhase = cellPhaseActiveIdxs(idx);
       
        cosmeticCenteringShift = floor((maxPixLength-allCellPixelLengths(idx))/2);
                
        for pixIdx = 1:allCellPixelLengths(idx)
            
            littleOutputMatrix(idx,pixIdx+cosmeticCenteringShift,1) = normalizedAllmeanY{theframes(cellPhase)}{thecellnos(cellPhase)}(pixIdx);
            littleOutputMatrix(idx,pixIdx+cosmeticCenteringShift,2) = normalizedAllmeanY{theframes(cellPhase)}{thecellnos(cellPhase)}(pixIdx);
            littleOutputMatrix(idx,pixIdx+cosmeticCenteringShift,3) = normalizedAllmeanY{theframes(cellPhase)}{thecellnos(cellPhase)}(pixIdx);
        
        end
        
        h8=figure(8); hold on;
        theSignal = normalizedAllmeanY{theframes(cellPhase)}{thecellnos(cellPhase)};
        meanSignal = mean(theSignal);
        theNormalizedSignal = theSignal./meanSignal;
        xvector=[1:allCellPixelLengths(idx)]+cosmeticCenteringShift;
        plot(xvector,theNormalizedSignal,'Color',1-[1 1 1].*(idx/nrActiveCellPhaseFrames));
        [peaki,peaky]=peakfinder(theNormalizedSignal);
        peakx=xvector(peaki);
        plot(peakx,peaky,'or');
                        
        maxPeak=max(peaky);
        
        % save for later use
        ringTimePlotx(end+1)=idx;
        ringTimePloty(end+1)=maxPeak;                
        
    end
    
    % save ring probability densitiess for later
    multipleringTimePlotx{end+1}=ringTimePlotx;
    multipleringTimePloty{end+1}=ringTimePloty;    
    
    % plot this ring probability density
    h9=figure(9); clf;     
    plot(ringTimePlotx,ringTimePloty,'-','LineWidth',3,'Color',mycolor);
    xlabel('Cell phase progression [frames]');
    ylabel('Highest peak / mean fluor value');
    MW_makeplotlookbetter(15);
    saveas(h9, [theoutputdir 'schnitzringintensityprogression'  num2str(schnitzIdx) '.tif']); 
    
    h8=figure(8); 
    saveas(h8, [theoutputdir 'schnitzintensityprofileprogression'  num2str(schnitzIdx) '.tif']); 
    
    h7=figure(7); clf; imshow(littleOutputMatrix);
    saveas(h7, [theoutputdir 'schnitzkymo'  num2str(schnitzIdx) '.tif']); 
end

%%
h11=figure(11); hold on;    


for i=1:numel(multipleringTimePlotx)    
    
    plot(multipleringTimePlotx{i},multipleringTimePloty{i},'-','LineWidth',1,'Color',mycolor);
    
end;

xlabel('Cell phase progression [frames]');
    ylabel('Highest peak / mean fluor value');
    MW_makeplotlookbetter(15);

saveas(h11, [theoutputdir 'ringTimingAll.tif']); 

%% ring occurence vs. time

h10=figure(10); clf; hold on;    
xlabel('Cell phase progression [normalized]');
ylabel('Highest peak / mean fluor value');
MW_makeplotlookbetter(15);

multipleringTimePlotxNormalized=cell(1,numel(multipleringTimePlotx));

for i=1:numel(multipleringTimePlotx)    
    cellCycleMax = max(multipleringTimePlotx{i}-1);
    
    multipleringTimePlotxNormalized{i} = (multipleringTimePlotx{i}-1)./cellCycleMax;
    plot(multipleringTimePlotxNormalized{i},multipleringTimePloty{i},'-','LineWidth',1,'Color',mycolor);
    
end;

[meanValuesForBins, binCenters]=binnedaveraging(multipleringTimePloty,multipleringTimePlotxNormalized,bins)

%plot(binCenters,meanValuesForBins,'ok','MarkerFaceColor','k','MarkerSize',10);
errorbar(binCenters,meanValuesForBins,errorValuesForBins,'ok','MarkerFaceColor','k','MarkerSize',10);

xlim([0,1])

saveas(h10, [theoutputdir 'ringTimingNormalized.tif']); 

%% more sanity checks

%load (['G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos' num2str(datasetIndex) 'crop\data\pos' num2str(datasetIndex) 'crop-Schnitz.mat']);







