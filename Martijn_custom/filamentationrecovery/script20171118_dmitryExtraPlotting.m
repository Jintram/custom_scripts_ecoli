
numberOflengths=numel(allL);
lengthsPx=arrayfun(@(i) numel(allD{i}), 1:numberOflengths);
lengthLongestCellPx=max(lengthsPx);

currentData=allD;

%%

% TODO: adjust to enable averaging

% Create prettier image array
prettyOutputImage = zeros(numberOflengths*2,max(lengthsPx));
% simpel output
for i = 1:numberOflengths
    
    %currentLength=max(allL{i});
    currentLength=allL(i);
    
    % normalize this data series
    cellProfile = currentData{i};
    cellProfile = cellProfile-min(cellProfile);
    % resize data to longest cell
    cellProfileResized = imresize(cellProfile,[1,lengthLongestCellPx],'Method','bilinear');

    % create output image
    % 
    prettyOutputImage(round((currentLength*2)-1),:) = cellProfileResized;
    prettyOutputImage(round((currentLength*2)),:) = cellProfileResized;
    
    %prettyOutputImage(i,:) = [cellProfile zeros(1,lengthLongestCell-numel(myData{i}))];
    %outputImage(i,:) = [cellProfile zeros(1,lengthLongestCell-numel(myData{i,:}))];
end

%output.(dataNames{dataSetIdx}).simpleOutputImage=simpleOutputImage;
%output.(dataNames{dataSetIdx}).prettyOutputImage=prettyOutputImage;


%% 
% copied from mw_getstatisticsandmakefigure

h2=figure(2); clf; 
imagesc(imrotate(prettyOutputImage,90)); 
hold on;

greenColorMap = makeColorMap([1 1 1],[65 148 68]./255);%,[105 189 69]./255)
redColorMap = makeColorMap([1 1 1],[230 30 37]./255);%,[230 30 37]./255)
colormap(greenColorMap);
%colorbar;

inputSettings.rangeIn = [size(output.D.prettyOutputImage,2),1];
inputSettings.desiredSpacing = .25;
inputSettings.rangeOut = [0,1];
[tickLocationsOldNewMetric, correspdongingLabels] = labelremapping(inputSettings);

set(gca,'XTick',[]);
set(gca,'YTick',tickLocationsOldNewMetric,'YTickLabel',correspdongingLabels);

xlabel('Length of cell [a.u.]');
ylabel(['Relative location along cell']);


%%

allD=allD;
allLPx = arrayfun(@(i) numel(allD{i}), 1:numel(allD));
allx = arrayfun(@(i) allL(i).*[1:allLPx(i)]./allLPx(i), 1:numel(allLPx),'UniformOutput',0);

counter=0;
allBacterialLengths=[];
%allFrameNrs=[6:13,14,16:22,24:32];
allxvalues          = {};
allyvalues          = {};
sumy                = [];
for idx = 1:numel(allL)


    counter=counter+1;

    WINDOWBORDERS   = [3:6:36];
    BARCOLORS       = [240,155,34; 45 177 65; 37 156 190; 131 84 162; 241 88 58]./255;

    %MW_straightenbacteria(p, [6:13,14,16:22,24:32], 'y') 


    distanceBacterium = allx{idx};
    bacterialLength = max(distanceBacterium);

    figure(2+counter); clf; hold on;
    plot(distanceBacterium, allD{idx},'o')

    myMaxY=max(allD{idx});

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
    allyvalues{end+1}          = allD{idx};
    sumy(end+1)                = sum(allD{idx});
    
    
end

%%
% Code copied from filamentationrecovery\script20170925_newMinProfiles3

figure(5); clf; hold on;
selectOnesWithPatternYesNo=0;

% make a color map for the signal strengths
myColorMapIntensities = makeColorMap([1 0 0],[0 1 0],[0 0 1],100);
sumyToMap = int32(99.*sumy./max(sumy)+1);

% init some parameters
[sortedLengths,sortedIndices] = sort(allBacterialLengths);

numelBacs = numel(sortedLengths);

sqrtNumel1 = ceil(sqrt(numelBacs));
sqrtNumel2 = sqrtNumel1;
if (sqrtNumel1.*(sqrtNumel1-1))>=numelBacs % trick to get 'm fit a bit tighter 
    sqrtNumel2=sqrtNumel1-1;
end

counter=0;
for plotIdx1 = 1:sqrtNumel1
    for plotIdx2 = 1:sqrtNumel2

        counter=counter+1;

        plotNr =  (plotIdx1-1)*sqrtNumel2+plotIdx2
        if plotNr>numelBacs
            break
        end

        % bacterium number (w. respect to defined by allFrameNrs)
        theBac = sortedIndices(plotNr);

        % color for signal of minD-YFP
        bacColor = myColorMapIntensities(sumyToMap(theBac),:);

        % plot the signal
        ax=subplot(sqrtNumel1,sqrtNumel2,plotNr);
        if ~selectOnesWithPatternYesNo                
            l=plot(allxvalues{theBac},allyvalues{theBac},'-','Color',bacColor);
            xlim([0,allBacterialLengths(theBac)]);
            title(['#' num2str(counter) ', L=' num2str(allBacterialLengths(theBac))]);% ', N=' num2str(regime)]);
        else
            % make a pretty plot
            l=plot(allxvalues{theBac},allyvalues{theBac}./max(allyvalues{theBac}),'-','Color',LINECOLOR,'LineWidth',2);
            xticks([0,round(max(allxvalues{theBac}))]);
            xlim([0,round(max(allxvalues{theBac}))]);
            ylim([0,1]);
            MW_makeplotlookbetter(20);
        end

        %regime=sum(WINDOWBORDERS<allBacterialLengths(theBac))+1;

        % cosmetics                        
        set(ax,'YTick',[]);

    end
end
