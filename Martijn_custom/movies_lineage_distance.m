
disp('Perhaps clear all?'); pause(3);

% Settings
% ===
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRange=[39:265]; % 39:265 for 2015-05-01/pos1crop 
myRange=[1:815]; % 1:815 for 2014_06_18/pos4crop 
% 39-200 is in 2015-05-01/pos8crop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputDir = 'U:\Matlab\Saved-analyses\';

MAKEIMAGES=1;

SMALLESTCELL = 0; % y-axis for cell length
BIGGESTCELL  = 6;

MAXGROWTHRATE = 2; % (dbls/hr)
ANCESTORDILUTION = .75; % How much color mixing dies out from subsequent ancestry


% Load schnitz stuff
%theName = ['pos1crop','2014-05-01']; p = DJK_initschnitz('pos1crop','2014-05-01','e.coli.AMOLF','rootDir','F:\A_Tans1_step4a_partially_analyzed_analysis\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
theName = ['pos4crop','2014_06_18']; p = DJK_initschnitz('pos4crop','2014_06_18','e.coli.AMOLF','rootDir','F:\A_Tans1_step4a_partially_analyzed_analysis\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
%theName = ['pos2crop','2014_06_18']; p = DJK_initschnitz('pos2crop','2014_06_18','e.coli.AMOLF','rootDir','F:\A_Tans1_step4a_partially_analyzed_analysis\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');

% Load the schnitzcells
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

% Now obtain the schnitzes that are living at the end of the movie
lastFrame = myRange(end);
[youngSchnitzes,cellTimePoints,cellLengths,timemin,timemax,lengthmin,lengthmax] =  ...
        givemeschnitzinfoforframe(p, schnitzcells, lastFrame);

% Relating all schnitzes to their final offspring can best be done going 
% backwards into the lineage.
YoungestOffspringForSchnitz = cell(size(schnitzcells)); % This list, with 
                % index i, contains the youngest schnitzes (i.e. end of 
                % lineage) that are grand-grand-..-children of schnitz i.
for youngGun = youngSchnitzes 
    % go back into ancestry until we have already seen this ancestor in
    % another search
    currentSchnitz = youngGun;
    done = 0;
    while ~done       
        if (isempty(currentSchnitz) ||  (currentSchnitz == 0))
            done = 1; % we're at the oldest schnitz in this lineage; done
        else
            % record that current schnitz has "youngGun" as final offspring
            YoungestOffspringForSchnitz{currentSchnitz}(end+1) = youngGun;
            % update current schnitz to its parent
            currentSchnitz = schnitzcells(currentSchnitz).P;
        end                
    end
end

% Now we make a colormap for each schnitz, comprising of a mix color of its
% offspring.
% First pick some colours for the young guns:
nrYoungGuns = numel(youngSchnitzes);
YoungGunColors = distinguishable_colors(nrYoungGuns, [0,0,0]);
% Now label each schnitz for its relatedness to the final offspring
IncestColor = zeros(size(schnitzcells,2),3);
for i = 1:numel(YoungestOffspringForSchnitz)   
    % only create a color for it when a youngest offspring has been
    % determined
    if isempty(YoungestOffspringForSchnitz{i}) break; end;

    % youngSchnitzes is a list of schnitznumbers, this search locates
    % the indices of that list for all schnitzes that are in 
    % YoungestOffspringForSchnitz{i}.
    colorIdx = find(ismember(youngSchnitzes,YoungestOffspringForSchnitz{i}));
        
    % Get all colors of the ancestry
    allYoungColors = YoungGunColors(colorIdx,:);
    if size(allYoungColors,1) > 1        
        IncestColor(i,:) = mean(allYoungColors);
    else
        IncestColor(i,:) = allYoungColors;
    end;
end

% When you don't want to mix children in:
% aColorMap = distinguishable_colors(numel(schnitzcells), [0,0,0]);


%% loop over the frames and obtain information (preprocessing)
%===
clear DATA; DATA(max(myRange)).a = 0;
for fr = myRange
    
    % Plot the image
    % ===

    % load the seg file (contains segmented and phase img)
    name= [p.segmentationDir,p.movieName,'seg',str3(fr)];
    load(name);        

    % Save Lc, note that p always remains the same
    DATA(fr).Lc = Lc;
    DATA(fr).phsub = phsub;
    
    % Get some info on this frame
    % ===
    [theSchnitzes,cellTimePoints,cellLengths,timemin,timemax,lengthmin,lengthmax,mapIndexes, output] = ...
      givemeschnitzinfoforframe(p, schnitzcells, fr);

    % Store this info
    DATA(fr).theSchnitzes   = theSchnitzes;
    DATA(fr).cellTimePoints = cellTimePoints;
    DATA(fr).cellLengths    = cellLengths;
    DATA(fr).mapIndexes     = mapIndexes;
    DATA(fr).output         = output;
    
    disp(['Analyzing range ' num2str(min(myRange)) '-' num2str(max(myRange)) ', now at frame '  num2str(fr)]);
end    

msgbox('Done');

%% Loop again, coloring them and creating plots
%==
if MAKEIMAGES

% figure
h=figure(1); clf; hold on; set(h, 'Position', [0,0,1200,500]);
h=figure(2); clf; hold on; set(h, 'Position', [0,0,1200,500]);
currentDateTime = datestr(now);
currentDateTime = strrep(currentDateTime, ':', '-');
currentDateTime = strrep(currentDateTime, ' ', '_');
for fr = myRange
    % Now we edit color map such that sisters look like their mothers, and a
    % bit less like their grandmothers, etc.
    % TODO: this code actually performs the same operation again and again;
    % this can be prevented..
    ancestryMap = ones(size(schnitzcells,2)+1,3);
    for i = DATA(fr).theSchnitzes
      colorsAncestry = []; currentAncestor = i; done = 0;
      relatedness = 1;
      while ~done
          colorsAncestry(end+1,:) = IncestColor(currentAncestor,:).*relatedness;
          currentAncestor = schnitzcells(currentAncestor).P;
          if (isempty(currentAncestor) || (currentAncestor == 0))
              done = 1;
          end
          relatedness = relatedness.*ANCESTORDILUTION;
      end
      if size(colorsAncestry,1) > 1
          ancestryMap(i,:) = sum(colorsAncestry)./sum(ANCESTORDILUTION.^[0:numel(colorsAncestry)-1]);        
      else
          ancestryMap(i,:) = colorsAncestry;
      end
    end

    % ABOVE BREAKS IT, BECAUSE IT IS DONE PER FRAME??

    % create a colormap for this frame
    % this means we need to find out for each shcnitz the number by which it
    % is known in this frame.  
    myCustomColorMap(1,:) = [0,0,0]; % note that way this works is that the 
                                   % value 0 is mapped to the first element
    for i = 1:numel(DATA(fr).mapIndexes)
      myCustomColorMap(DATA(fr).mapIndexes(i)+1,:) = ancestryMap(DATA(fr).theSchnitzes(i),:);
    end    

    % Let's also create a color map for growth rate
    growthRateColors = zeros(size(schnitzcells,2)+1,3);
    for i = 1:numel(DATA(fr).mapIndexes)
      colorintensity = DATA(fr).output{i}.muP11_all/MAXGROWTHRATE;
      growthRateColors(DATA(fr).mapIndexes(i)+1,:) = [colorintensity, 0, 0]; % +1 for black at start
    end    

    % make image of cells 
    p.showPerim = 0;
    outim = PN_imshowlabel(p, DATA(fr).Lc, 0,0,0, 'customColors', myCustomColorMap);%,'phaseImage',); % note i'm feeding the same img as "previous" img
    % crop img
    %outim = imcrop(outim, [395   225   650   650]);

    % growthrates image
    p.showPerim = 0;
    imGrowthRates = PN_imshowlabel(p, DATA(fr).Lc, 0,0,0, 'customColors', growthRateColors);%,'phaseImage',); % note i'm feeding the same img as "previous" img

    % Outline image
    p.showPerim = 1;
    imOutlines = PN_imshowlabel(p, DATA(fr).Lc, 0,0,0,'customColors', myCustomColorMap, 'phaseImage', DATA(fr).phsub); % note i'm feeding the same img as "previous" img

    % Plot image of cells and length plot
    % ===

    change_current_figure(1); h=1;

    %  image of cells
    subplottight(1,3,3); 
    imshow(outim,[]);
    text(20,20,['frame ' sprintf('%05d', fr)],'Color','k','FontWeight','bold','BackgroundColor','white');  

    %  image of cells
    subplottight(1,3,2); 
    imshow(imGrowthRates,[]);
    text(20,20,['frame ' sprintf('%05d', fr)],'Color','k','FontWeight','bold','BackgroundColor','white');  

    %  image of cells
    subplottight(1,3,1); 
    imshow(imOutlines,[]);
    text(20,20,['frame ' sprintf('%05d', fr)],'Color','k','FontWeight','bold','BackgroundColor','white');  
    
    theDir = [outputDir 'movie-' theName '-' currentDateTime '\'];
    if ~exist(theDir)
        mkdir(theDir)
    end
    saveas(h, [theDir 'movietest_cells_' sprintf('%05d',fr) '.jpg']); % below is another save command for the other fig

    % their lengths in a plot
    change_current_figure(2); h = 2;
    for j = 1:numel(DATA(fr).theSchnitzes)
    plot(fr, DATA(fr).cellLengths(j), '.','Color',myCustomColorMap(DATA(fr).mapIndexes(j)+1,:)) % +1 because 1st color = black    
    hold on;
    %plot(fr, cellLengths(j), '.','Color',ind2rgb(j,colormap(lines)))
    end
    text(20,20,['frame ' num2str(fr)],'Color','k','FontWeight','bold','BackgroundColor','white');
    xlim([min(myRange),max(myRange)]), ylim([SMALLESTCELL,BIGGESTCELL]); % redundant  

    title('Cell lengths','FontSize',15)
    xlabel('Frame number','FontSize',15)
    ylabel('Cell length ({\mu}m)','FontSize',15)    

    saveas(h, [theDir 'movietest_lengthts_' sprintf('%05d',fr) '.jpg']);
    
    disp(['Analyzing range ' num2str(min(myRange)) '-' num2str(max(myRange)) ', now at frame '  num2str(fr)]);
    
end

msgbox('Done');

end

%% Determine correlations between all cell pairs in growth rate
% ===
NEIGHBORDISTANCE = 40; % not sure which units this is, prob. pixels
RELATEDNESSCUTOFF = 1; % 1 = sisters, 2 = same granny, etc

correlationsPerFrame = []; allpairsN = []; 
correlationsPerFrameNeighbors = [];neighborN = [];
correlationsPerFrameSisters = []; sisterN = [];
correlationsPerFrameUnrelatedNeighbors = []; unrelatedNeighborN = [];
for fr = myRange
    
    allrelations = []; % TODO remove
    
    theGrowthRates = [];
    for i=1:numel(DATA(fr).output)
        theGrowthRates(end+1) = DATA(fr).output{i}.muP11_fitNew_all;
    end % TODO this should be possible in one command?
    
    % theoretically, all pairs should not have a correlation in growth rate
    theGrowthRates=theGrowthRates./mean(theGrowthRates); % normalize
    theGrowthRates=theGrowthRates-mean(theGrowthRates);
    N=numel(theGrowthRates); nrPairs = N*(N-1)/2; n = 0;
    pairsOfAll = zeros(nrPairs,2);
    pairsOfNeighbors = []; pairsOfRelated = []; pairsOfUnrelatedNeighbors = [];
    for i = 1:(N-1)
        for j = (i+1):N %i:N

            if i~=j        
                n = n + 1;
                % All pairs
                pairsOfAll(n,:) = [theGrowthRates(i),theGrowthRates(j)];
                
                % Create distance criterium subset
                % obtain distance between bacteria
                distance = sqrt( ...
                    (DATA(fr).output{i}.cenx_cent - DATA(fr).output{j}.cenx_cent)^2+ ...
                    (DATA(fr).output{i}.ceny_cent - DATA(fr).output{j}.ceny_cent)^2  ...
                                );
                       
                % create pairs of neighbors
                if distance<NEIGHBORDISTANCE                    
                    pairsOfNeighbors = [pairsOfNeighbors; theGrowthRates(i),theGrowthRates(j)];
                    % disp(['hit in fr ' num2str(fr)]);
                end
                
                % create pairs of sisters, note self-interaction is
                % excluded above
                relatedness = distancelca(p, schnitzcells, DATA(fr).theSchnitzes(i), DATA(fr).theSchnitzes(j), RELATEDNESSCUTOFF+1);                
                allrelations(end+1) = relatedness; % TODO remove
                if relatedness<=RELATEDNESSCUTOFF 
                    pairsOfRelated = [pairsOfRelated; theGrowthRates(i),theGrowthRates(j)];
                end
                
                % neighbors but not related
                if (distance<NEIGHBORDISTANCE) && (relatedness>RELATEDNESSCUTOFF || isnan(relatedness))
                    pairsOfUnrelatedNeighbors = [pairsOfUnrelatedNeighbors; theGrowthRates(i),theGrowthRates(j)];
                end
                
                
                
            end

        end    
    end

    % Store correlations for all pairs
    correlationsPerFrame(end+1) = corr(pairsOfAll(:,1),pairsOfAll(:,2),'type','Pearson');
    
    % Calculate and store correlations for neighbouring cells
    if size(pairsOfNeighbors,1) == 0
        correlationsPerFrameNeighbors(end+1) = NaN;
    else
        correlationsPerFrameNeighbors(end+1) = corr(pairsOfNeighbors(:,1),pairsOfNeighbors(:,2),'type','Pearson');
    end
    
    % Calculate and store correlations for sister cells
    if size(pairsOfRelated,1) == 0
        correlationsPerFrameSisters(end+1) = NaN;
    else
        correlationsPerFrameSisters(end+1) = corr(pairsOfRelated(:,1),pairsOfRelated(:,2),'type','Pearson');
    end
    
    
    % Calculate and store correlations for sister cells
    if size(pairsOfUnrelatedNeighbors,1) == 0
        correlationsPerFrameUnrelatedNeighbors(end+1) = NaN;
    else
        correlationsPerFrameUnrelatedNeighbors(end+1) = corr(pairsOfUnrelatedNeighbors(:,1),pairsOfUnrelatedNeighbors(:,2),'type','Pearson');
    end
    
    % counting the pairs
    allpairsNforFrame = size(pairsOfAll,1); 
    allpairsN(end+1) = allpairsNforFrame;
    
    neighborNforFrame = size(pairsOfNeighbors,1);
    neighborN(end+1) = neighborNforFrame;
    
    sisterNforFrame = size(pairsOfRelated,1);
    sisterN(end+1) = sisterNforFrame;
    
    unrelatedNeighborsNforFrame = size(pairsOfUnrelatedNeighbors,1);
    unrelatedNeighborN(end+1) = unrelatedNeighborsNforFrame;
    
    disp(['Selected ' num2str(neighborNforFrame) '/' num2str(allpairsNforFrame) ' as neighbors']);
    disp(['Selected ' num2str(sisterNforFrame) '/' num2str(allpairsNforFrame) ' as related']);
    disp(['Selected ' num2str(unrelatedNeighborsNforFrame) '/' num2str(allpairsNforFrame) ' as unrelated neighbors']);
    
    
    % Now select pairs that have certain distance
    
    
    %{
    std1 = std(pairsOfLengths(:,1));
    std2 = std(pairsOfLengths(:,2));
    mycorrManual = mean((pairsOfLengths(:,1).*pairsOfLengths(:,2)))/(std1*std2)
    %}

    disp(['Analyzing range ' num2str(min(myRange)) '-' num2str(max(myRange)) ', now at frame '  num2str(fr)]);
    
end

msgbox('Done');

% Save data
theDir = [outputDir 'output\'];
currentDateTime = datestr(now);
currentDateTime = strrep(currentDateTime, ':', '-');
currentDateTime = strrep(currentDateTime, ' ', '_');
if ~exist(theDir)
    mkdir(theDir)
end
save( [theDir theName 'calculationoutput-correlationsneighbors_' currentDateTime '.mat'] );
disp(['Data saved to ' theDir]);

%% plot
DESIREDNPAIRSFORSTATS = 10;

figure (3); clf;
subplot(1,3,1); hold on;
plot([min(myRange),max(myRange)],[0,0],'k');
l1=plot(myRange,correlationsPerFrame,'k','Linewidth',2);
l2=plot(myRange,correlationsPerFrameNeighbors,'b','Linewidth',2);
l3=plot(myRange,correlationsPerFrameSisters,'r','Linewidth',2);
l4=plot(myRange,correlationsPerFrameUnrelatedNeighbors,'g','Linewidth',2);
ylim([-1,1]);
xlim([min(myRange),max(myRange)])
xlabel('Frame #'); ylabel('Correlation');
%legend(l1,'All pairs');

subplot(1,3,2); hold on;
plot([min(neighborN),max(neighborN)],[0,0],'-k');
l1=plot(neighborN,correlationsPerFrameNeighbors,'bs','Linewidth',2);
l2=plot(sisterN,correlationsPerFrameSisters,'r^','Linewidth',2);
l3=plot(unrelatedNeighborN,correlationsPerFrameUnrelatedNeighbors,'vg','Linewidth',2);
l4=plot(allpairsN,correlationsPerFrame,'ko','Linewidth',4);
ylim([-1,1]);
xlabel('Number of pairs considered'); ylabel('Correlation');
xlim([min(neighborN),max(neighborN)])
legend([l1,l2,l3,l4],{'neighbors', 'sisters', 'unrelated neighbors', 'all'})

subplot(1,3,3); hold on;
l1=plot(myRange,allpairsN,'k--','Linewidth',3);
l2=plot(myRange,neighborN,'b-','Linewidth',3);
l3=plot(myRange,sisterN,'r-','Linewidth',3);
l4=plot(myRange,unrelatedNeighborN,'g-','Linewidth',3);
xlabel('Number of pairs considered'); ylabel('Correlation');
xlim([min(myRange),max(myRange)])
xlabel('Frame #'); ylabel('Number of pairs considered');
%legend(l1,'All pairs');

% Calculate mean lines; for neighbors
theMeansNeighbors = []; theStdsNeighbors = [];
uniqueSortedNeighborN = sort(unique(neighborN));
for i = uniqueSortedNeighborN
    theMeansNeighbors(end+1) = mean(correlationsPerFrameNeighbors(neighborN(:)==i))
    theStdsNeighbors(end+1) = std(correlationsPerFrameNeighbors(neighborN(:)==i))
end
% Calculate mean lines; for sisters
theMeansSisters = []; theStdsSisters = [];
uniqueSortedsisterN = sort(unique(sisterN));
for i = uniqueSortedsisterN
    theMeansSisters(end+1) = mean(correlationsPerFrameSisters(sisterN(:)==i))
    theStdsSisters(end+1) = std(correlationsPerFrameSisters(sisterN(:)==i))
end
% Calculate mean lines; for unrelated neighbors
theMeansUnrelatedNeighbors = []; theStdsUnrelatedNeighbors = [];
uniqueSortedUnrelatedNeighborN = sort(unique(unrelatedNeighborN));
for i = uniqueSortedUnrelatedNeighborN
    theMeansUnrelatedNeighbors(end+1) = mean(correlationsPerFrameUnrelatedNeighbors(unrelatedNeighborN(:)==i))
    theStdsUnrelatedNeighbors(end+1) = std(correlationsPerFrameUnrelatedNeighbors(unrelatedNeighborN(:)==i))
end
% Calculate mean lines; for all pairs
theMeansAllPairs = []; theStdsAllPairs = [];
uniqueSortedallpairsN = sort(unique(allpairsN));
for i = uniqueSortedallpairsN
    theMeansAllPairs(end+1) = mean(correlationsPerFrame(allpairsN(:)==i))
    theStdsAllPairs(end+1) = std(correlationsPerFrame(allpairsN(:)==i))
end

figure(4), clf,  hold on;
xlim([min([uniqueSortedNeighborN, uniqueSortedsisterN,uniqueSortedUnrelatedNeighborN ]),max([uniqueSortedNeighborN, uniqueSortedsisterN,uniqueSortedUnrelatedNeighborN ])]);
ylim([-1,1])
l1=errorbar(uniqueSortedNeighborN,theMeansNeighbors,theStdsNeighbors,'ob','Linewidth',2)
errorbar_tick(l1)
l2=errorbar(uniqueSortedsisterN,theMeansSisters,theStdsSisters,'or','Linewidth',2)
errorbar_tick(l2)
l3=errorbar(uniqueSortedUnrelatedNeighborN,theMeansUnrelatedNeighbors,theStdsUnrelatedNeighbors,'og','Linewidth',2)
errorbar_tick(l3)
l4=errorbar(uniqueSortedallpairsN,theMeansAllPairs,theStdsAllPairs,'ok','Linewidth',2)
errorbar_tick(l4)
legend([l1,l2,l3,l4],{'neighbors','sisters','unrelated neighbors','all pairs'},'Location','southeast')
plot([min([uniqueSortedNeighborN, uniqueSortedsisterN,uniqueSortedUnrelatedNeighborN ]),max([uniqueSortedNeighborN, uniqueSortedsisterN,uniqueSortedUnrelatedNeighborN ])],[0,0],'-k')
title('Correlations between pairs')
xlabel('Number of pairs considered')
ylabel('Correlation (normalized)')
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal')
set(gca,'FontSize',15)

figure(5); clf; hold on;
xlim([min([uniqueSortedallpairsN ]),max([uniqueSortedallpairsN ])]);
ylim([-1,1])
l4=errorbar(uniqueSortedallpairsN,theMeansAllPairs,theStdsAllPairs,'ok','Linewidth',2)
errorbar_tick(l4)
plot([min([uniqueSortedallpairsN ]),max([uniqueSortedallpairsN ])],[0,0],'-k')
title('Correlations between all pairs')
xlabel('Number of pairs considered')
ylabel('Correlation (normalized)')
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal')
set(gca,'FontSize',15)

% Mean values
% Calculate mean lines; for neighbors
theMeanNeighbors = mean(correlationsPerFrameNeighbors(neighborN(:)>DESIREDNPAIRSFORSTATS))
theStdNeighbors = std(correlationsPerFrameNeighbors(neighborN(:)>DESIREDNPAIRSFORSTATS))

% Calculate mean lines; for sisters
theMeanSisters = mean(correlationsPerFrameSisters(sisterN(:)>DESIREDNPAIRSFORSTATS))
theStdSisters = std(correlationsPerFrameSisters(sisterN(:)>DESIREDNPAIRSFORSTATS))

% Calculate mean lines; for unrelated neighbors
theMeanUnrelatedNeighbors = mean(correlationsPerFrameUnrelatedNeighbors(unrelatedNeighborN(:)>DESIREDNPAIRSFORSTATS))
theStdUnrelatedNeighbors = std(correlationsPerFrameUnrelatedNeighbors(unrelatedNeighborN(:)>DESIREDNPAIRSFORSTATS))

% Calculate mean lines; for all pairs
theMeanAllPairs = mean(correlationsPerFrame(allpairsN(:)>DESIREDNPAIRSFORSTATS))
theStdAllPairs = std(correlationsPerFrame(allpairsN(:)>DESIREDNPAIRSFORSTATS))

figure(6); clf;
barwitherr([theStdNeighbors,theStdSisters,theStdUnrelatedNeighbors,theStdAllPairs],[theMeanNeighbors,theMeanSisters,theMeanUnrelatedNeighbors,theMeanAllPairs])
Labels = {'Neighbors', 'Related', 'Unr. neighbors', 'All pairs'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal')
set(gca,'FontSize',15)
ylabel('Correlation (normalized)')
title('Correlations between cell pairs in a colony')

%% 
myFrame = 100;
xdata=[]; ydata=[];
for i = 1:numel(DATA(myFrame).output)
    i
    xdata(end+1) = DATA(myFrame).output{i}.cenx_cent;
    ydata(end+1) = DATA(myFrame).output{i}.ceny_cent;    
end
figure, plot(xdata, ydata,'x');

%{
currentFolder = pwd;
outputDir = [currentFolder '\output\'];
currentDateTime = datestr(now);
currentDateTime = strrep(currentDateTime, ':', '-');
currentDateTime = strrep(currentDateTime, ' ', '_');
if ~exist(outputDir)
    mkdir(outputDir)
end
save( [outputDir 'mydata_' currentDateTime '.mat'] );

disp( ['Data saved to ' outputDir]);
%}

%% 




%{
figure(1), imshow(Lc,[])
figure(2), imshow(phsub,[])

figure(2), imshow(outim,[])
%}
