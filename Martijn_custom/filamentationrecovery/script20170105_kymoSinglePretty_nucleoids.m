
%OUTPUTDIR = 'H:\EXPERIMENTAL_DATA_2017\2017-10-12_hupA-mRuby2\pos1cropa2\analysis\straightenedCells\kymographsinglecell\';
OUTPUTDIR = [p.analysisDir 'straightenedCells\kymographsinglecell\'];

if ~exist('OUTPUTDIR','dir')
    mkdir(OUTPUTDIR);
end

OUTLINECOLOR = [48, 197, 221]/255;

%% Load data

%load('G:\EXPERIMENTAL_DATA_2016\2016-04-07_asc777_temperatureRecovery\pos2crop\analysis\straightenedCells\2016-04-07pos2crop_straightFluorData.mat')
%load('G:\EXPERIMENTAL_DATA_2016\2016-04-07_asc777_temperatureRecovery\pos2crop\data\pos2crop-Schnitz.mat');
%load('H:\EXPERIMENTAL_DATA_2017\2017-10-12_hupA-mRuby2\pos1cropa2\data\pos1cropa2-Schnitz.mat');
load([p.tracksDir p.movieName '-Schnitz.mat']);

%load('H:\EXPERIMENTAL_DATA_2017\2017-10-12_hupA-mRuby2\pos1cropa2\analysis\straightenedCells\2017-10-12pos1cropa2_straightFluorData.mat');
load([p.analysisDir 'straightenedCells\' p.movieDate p.movieName '_straightFluorData.mat']);

micronsPerPixel = 0.0431; % see G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos2crop_parameters.mat, p struct

%% backwards compatibility
if ~exist('allmeanY','var')
    allmeanY=allrawFluor;
    warning('allmeanY created from allrawFluor');
end

%% normalization of signals

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

%% Create kymograph of desired schnitz
%{
if ~exist('SCHNITZIDX','var')
    SCHNITZIDX=11; % SCHNITZIDX=11;
end
%}

for SCHNITZIDX = 1:numel(schnitzcells)

    %%
    SCHNITZIDX
    
    % for ring timing plots
    ringTimePlotx=[];    
    ringTimePloty=[];

    % frame and cellnos applicable for this schnitz needed for calculations
    theframes  = schnitzcells(SCHNITZIDX).frame_nrs;
    thecellnos = schnitzcells(SCHNITZIDX).cellno;

    cellLifeLength = numel(theframes); % in frames, for administration

    allCellPixelLengths=[]; cellPhaseActiveIdxs = [];
    for i = 1:numel(theframes)
        try
            if ~isempty(allmeanY{theframes(i)})             
                    allCellPixelLengths(end+1) = numel(allmeanY{theframes(i)}{thecellnos(i)});
                    cellPhaseActiveIdxs(end+1) = i;
            end
        catch
            warning('Warning, something went wrong..');            
        end
    end

    maxPixLength = max(allCellPixelLengths);
    nrActiveCellPhaseFrames = numel(cellPhaseActiveIdxs);

    if nrActiveCellPhaseFrames==0
        disp(['Skipping schnitz ' num2str(SCHNITZIDX) ' since nrActiveCellPhaseFrames is 0.']);
        continue
    end
    
    littleOutputMatrix = ones(nrActiveCellPhaseFrames,maxPixLength,3); %% 
    %littleOutputMatrix(:,:,3) = ones(nrActiveCellPhaseFrames,maxPixLength);

    bacterialOutline1=zeros(1,nrActiveCellPhaseFrames*2); bacterialOutline2=zeros(1,nrActiveCellPhaseFrames*2); xvalues=zeros(1,nrActiveCellPhaseFrames*2);
    for idx = 1:nrActiveCellPhaseFrames        

        cellPhase = cellPhaseActiveIdxs(idx);

        cosmeticCenteringShift = floor((maxPixLength-allCellPixelLengths(idx))/2);

        for pixIdx = 1:allCellPixelLengths(idx)

            littleOutputMatrix(idx,pixIdx+cosmeticCenteringShift,1) = normalizedAllmeanY{theframes(cellPhase)}{thecellnos(cellPhase)}(pixIdx);
            littleOutputMatrix(idx,pixIdx+cosmeticCenteringShift,2) = normalizedAllmeanY{theframes(cellPhase)}{thecellnos(cellPhase)}(pixIdx);
            littleOutputMatrix(idx,pixIdx+cosmeticCenteringShift,3) = normalizedAllmeanY{theframes(cellPhase)}{thecellnos(cellPhase)}(pixIdx);

        end

        % Saving where the edges are to color them later
        % The elaborate indexing is to really outline the pixels, ie. their
        % edges
        xvalues(2*idx-1)=idx-.5; % left edge of pixel
        xvalues(2*idx)  =idx+.5;   % right edge of pixel
        pixIdx=1;
        bacterialOutline1(2*idx-1) = pixIdx+cosmeticCenteringShift;    
        bacterialOutline1(2*idx)   = pixIdx+cosmeticCenteringShift;    
        pixIdx=allCellPixelLengths(idx);
        bacterialOutline2(2*idx-1) = pixIdx+cosmeticCenteringShift+1;
        bacterialOutline2(2*idx)   = pixIdx+cosmeticCenteringShift+1;

    end

    %%
    h=figure(1); clf; hold on;
    NRXTICKS=10; NRYTICKS=10;

    littleOutputMatrixRotated = imrotate(littleOutputMatrix,90);

    ax=image(imrotate(littleOutputMatrix,90));
    %ax=image(littleOutputMatrixRotated);
    %ax=image(littleOutputMatrix);

    % Calculate x labels
    TIMEDISTANCE=60;
    maxTime = max(schnitzcells(SCHNITZIDX).time)-min(schnitzcells(SCHNITZIDX).time);
    maxPixelsTime = size(littleOutputMatrix,1);
    timePerPixel = maxPixelsTime/maxTime;
    xticksInHrs    = [0:TIMEDISTANCE:maxTime];
    xticksInPixels = xticksInHrs.*timePerPixel;
    %rangeElementsX=linspace(1,nrActiveCellPhaseFrames,NRXTICKS);
    %rangeTimepoints=linspace(0,cellLifeTime,NRXTICKS);

    myxlabels = {};
    for i=1:numel(xticksInHrs)
        myxlabels{i}= sprintf('%.0f', xticksInHrs(i));
    end

    % Calculate Y labels
    ABSOLUTEMICRONDISTANCE=10;
    distanceYticksPixels = ABSOLUTEMICRONDISTANCE/micronsPerPixel;
    %rangeElementsY=linspace(1,size(littleOutputMatrix,2),NRYTICKS);
    maxValuePixels=size(littleOutputMatrix,2);
    maxValueMicrons=maxValuePixels*micronsPerPixel;
    %rangeYpoints=linspace(cellMaxSize,0,NRYTICKS);%-.5*cellMaxSize;
    yticksInMicrons=[maxValuePixels:-distanceYticksPixels:0];
    yticksInPixels=[0:ABSOLUTEMICRONDISTANCE:maxValueMicrons];

    myylabels = {};
    for i=1:numel(yticksInMicrons)
        myylabels{i}= sprintf('%.0f', yticksInPixels(i));
    end

    set(ax.Parent,'XTick',xticksInPixels,'XTickLabel',myxlabels)
    set(ax.Parent,'YTick',yticksInMicrons(end:-1:1),'YTickLabel',myylabels(end:-1:1))
    %set(ax(1),'YTick',[0:0.2:1], 'XTickLabel',VALUES)

    % plot outlines
    plot(xvalues,bacterialOutline1,'-','Color',OUTLINECOLOR,'LineWidth',2);
    plot(xvalues,bacterialOutline2,'-','Color',OUTLINECOLOR,'LineWidth',2);

    xlim([.5, nrActiveCellPhaseFrames+.5]);

    xlabel('Time (min)');
    %ylabel(['FtsA-YFP signal' 10 'along cellular axis (a.u.)']);
    ylabel('Cellular axis (µm)')

    ylim([0,size(littleOutputMatrix,2)]) % not sure why needed, but addresses issue

    MW_makeplotlookbetter(20);

    %h.PaperUnits = 'centimeters';
    %h.PaperPosition = [0 0 13.1 6.6];

    if ~exist('NOSAVEPLEASE','var')   
        saveas(h,[OUTPUTDIR 'TIF_kymo_schnitz' num2str(SCHNITZIDX) '.tif'])
        saveas(h,[OUTPUTDIR 'SVG_kymo_schnitz' num2str(SCHNITZIDX) '.svg'])
        saveas(h,[OUTPUTDIR 'FIG_kymo_schnitz' num2str(SCHNITZIDX) '.fig'])
    end

    hFig3b=h;

end

disp('All done');


