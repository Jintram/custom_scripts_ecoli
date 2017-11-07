
%% Load data

load('H:\EXPERIMENTAL_DATA_2017\2017-09-19_MinD-YFP\PART4_MinD-YFP_movies\pos1\analysis\straightenedCells\PART4_MinDpos1_straightFluorData.mat');
load('H:\EXPERIMENTAL_DATA_2017\2017-09-19_MinD-YFP\PART4_MinD-YFP_movies\pos1\data\pos1-skeletonData.mat');

if ~exist('p','var') % if ~isfield(p,'micronsPerPixel')
    p.micronsPerPixel=   0.0438;
end


%% plotting all

allxvalues={}; allyvalues={}; sumy=[];

% convenient to look one up:
% PLOTNR=2; allFrameNrs(sortedIndices(PLOTNR))

counter=0;
allBacterialLengths=[];
%allFrameNrs=[6:13,14,16:22,24:32];
allFrameNrs=[6:10,12,13,16:22,24:32]; % filtered out bad ones: 11, 14
for frameNr = allFrameNrs


    counter=counter+1;

    WINDOWBORDERS   = [3:6:36];
    BARCOLORS       = [240,155,34; 45 177 65; 37 156 190; 131 84 162; 241 88 58]./255;

    %MW_straightenbacteria(p, [6:13,14,16:22,24:32], 'y') 


    distanceBacterium = allskeletonDistance{frameNr}{1}.*p.micronsPerPixel;
    bacterialLength = max(distanceBacterium);

    figure(2+counter); clf; hold on;
    plot(distanceBacterium, allrawFluor{frameNr}{1},'o')

    myMaxY=max(allrawFluor{frameNr}{1});

    if myMaxY>0
        ylim([0,myMaxY*1.1])
    end


    % this could be done prettier, but draws lines where expected
    if bacterialLength < WINDOWBORDERS(2)
        plot([.5 .5].*bacterialLength,[0,myMaxY*1.1],'-','Color',BARCOLORS(1,:))
    elseif bacterialLength < WINDOWBORDERS(3)
        for i=1:2:4
            x = i*[1/4 1/4].*bacterialLength
            plot(x,[0,myMaxY*1.1],'-','Color',BARCOLORS(2,:))
        end
    elseif bacterialLength < WINDOWBORDERS(4)
        for i=1:2:6
            x = i*[1/6 1/6].*bacterialLength
            plot(x,[0,myMaxY*1.1],'-','Color',BARCOLORS(3,:))
        end
    end

    % summary stats
    allBacterialLengths(end+1) = bacterialLength;
    allxvalues{end+1}          = distanceBacterium;
    allyvalues{end+1}          = allrawFluor{frameNr}{1};
    sumy(end+1)                = sum(allrawFluor{frameNr}{1});
    
    
end

%% Make combined plot

figure(101); clf; hold on;

% make a color map for the signal strengths
myColorMapIntensities = makeColorMap([1 0 0],[0 1 0],[0 0 1],100);
sumyToMap = int32(99.*sumy./max(sumy)+1);

% init some parameters
[sortedLengths,sortedIndices] = sort(allBacterialLengths);
 
numelBacs = numel(sortedLengths);

sqrtNumel = ceil(sqrt(numelBacs));

for plotIdx1 = 1:sqrtNumel
    for plotIdx2 = 1:sqrtNumel
    
        
        plotNr =  (plotIdx1-1)*sqrtNumel+plotIdx2
        if plotNr>numelBacs
            break
        end
        
        % bacterium number (w. respect to defined by allFrameNrs)
        theBac = sortedIndices(plotNr);
       
        % color for signal of minD-YFP
        bacColor = myColorMapIntensities(sumyToMap(theBac),:);
        
        % plot the signal
        ax=subplot(sqrtNumel,sqrtNumel,plotNr);
        l=plot(allxvalues{theBac},allyvalues{theBac},'-','Color',bacColor);
       
        %regime=sum(WINDOWBORDERS<allBacterialLengths(theBac))+1;
        
        % cosmetics
        title(['L=' num2str(allBacterialLengths(theBac))]);% ', N=' num2str(regime)]);
        xlim([0,allBacterialLengths(theBac)])
       
        set(ax,'YTick',[]);
        
    end
end

figure(102); colormap(myColorMapIntensities); colorbar;

%% 
relativeL=50; N=50;
pixFactor=2;

myMatrixWithResults = zeros(relativeL,N);

for frameNr = [6:13,14,16:22,24:32]

    %%    
    
    yvalues=allrawFluor{frameNr}{1};
    yvalues = (yvalues-min(yvalues))./(max(yvalues)-min(yvalues));
    
    distanceBacterium = allskeletonDistance{frameNr}{1}.*p.micronsPerPixel;
    bacterialLength = max(distanceBacterium);
    distanceBacteriumRelativeL = distanceBacterium./bacterialLength;
    
    [meanValuesForBins, binCenters]= binnedaveraging({distanceBacteriumRelativeL.*relativeL}, {yvalues}, 0:relativeL);
    
    %figure(2+counter); clf; hold on;
    %plot(binCenters, meanValuesForBins ,'o')
    
    %% put in matrix
    
    myMatrixWithResults(:,pixFactor.*round(bacterialLength)) = ...
        myMatrixWithResults(:,round(bacterialLength)) + meanValuesForBins';
    
end

% normalize matrix

%% make plot

h1=figure(100); clf;
%imagesc(imrotate(myMatrixWithResults,90)); 
imagesc(myMatrixWithResults); 

greenColorMap = makeColorMap([1 1 1],[65 148 68]./255);%,[105 189 69]./255)
redColorMap = makeColorMap([1 1 1],[230 30 37]./255);%,[230 30 37]./255)
colormap(redColorMap);
%colorbar;

%{
% recalculate y-axis
inputSettings.rangeIn = [size(myMatrixWithResults,2),1];
inputSettings.desiredSpacing = .25;
inputSettings.rangeOut = [0,1];
[tickLocationsOldNewMetric, correspdongingLabels] = labelremapping(inputSettings);

set(gca,'XTick',[]);
set(gca,'YTick',tickLocationsOldNewMetric,'YTickLabel',correspdongingLabels);

xlabel('Length of cell [a.u.]');
ylabel(['Relative location along cell']);
%}

