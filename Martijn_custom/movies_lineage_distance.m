
% Settings
% ===
MAXGROWTHRATE = 2; % (dbls/hr)
ANCESTORDILUTION = .75; % How much color mixing dies out from subsequent ancestry
myRange=[39:200] % 39-200 is in 2015-05-01/pos8crop 

SMALLESTCELL = 0; % y-axis for cell length
BIGGESTCELL  = 6;

% Load schnitz stuff
p = DJK_initschnitz('pos8crop','2014-05-01','e.coli.AMOLF','rootDir','F:\A_Tans1_step4a_partially_analyzed_analysis\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');

% For plotting lengths
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

%% Loop again, coloring them and creating plots
%==

% figure
h=figure(1), set(h, 'Position', [0,0,800,500]);
h=figure(2), set(h, 'Position', [0,0,800,500]);
h=figure(1); clf; hold on;
h=figure(2); clf; hold on;
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

    saveas(h, ['D:\Local_Playground\mymovietest\movietest_cells_' sprintf('%05d',fr) '.jpg']);

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

    saveas(h, ['D:\Local_Playground\mymovietest\movietest_lengths_' sprintf('%05d',fr) '.jpg']);
    
    disp(['Analyzing range ' num2str(min(myRange)) '-' num2str(max(myRange)) ', now at frame '  num2str(fr)]);
    
end

msgbox('Done');

%% Determine correlations between all cell pairs
% ===
NEIGHBORDISTANCE = 75; % not sure which units this is, prob. pixels

correlationsPerFrame = []; allpairsN = []; 
correlationsPerFrameNeighbors = [];neighborN = [];
for fr = myRange
        
    theLengths = DATA(fr).cellLengths;
    % theoretically, all length pairs should not have a correlation
    theLengths=theLengths./mean(theLengths); % normalize
    theLengths=theLengths-mean(theLengths);
    N=numel(theLengths); nrPairs = N*(N-1); n = 0;
    pairsOfLengths = zeros(nrPairs,2);
    pairsOfLengthsNeighbors = [];
    for i = 1:(N)
        for j = 1:N %i:N

            if i~=j        
                n = n + 1;
                % All pairs
                pairsOfLengths(n,:) = [theLengths(i),theLengths(j)];
                
                % Create distance criterium subset
                % obtain distance between bacteria
                distance = sqrt( ...
                    (DATA(fr).output{i}.cenx_cent - DATA(fr).output{j}.cenx_cent)^2+ ...
                    (DATA(fr).output{i}.ceny_cent - DATA(fr).output{j}.ceny_cent)^2  ...
                                );
                            
                if distance<NEIGHBORDISTANCE                    
                    pairsOfLengthsNeighbors = [pairsOfLengthsNeighbors; theLengths(i),theLengths(j)];
                    % disp(['hit in fr ' num2str(fr)]);
                end                            
                % 
                
                
            end

        end    
    end

    % Store correlations;
    correlationsPerFrame(end+1) = corr(pairsOfLengths(:,1),pairsOfLengths(:,2),'type','Pearson');
    if size(pairsOfLengthsNeighbors,1) == 0
        correlationsPerFrameNeighbors = NaN;
    else
        correlationsPerFrameNeighbors(end+1) = corr(pairsOfLengthsNeighbors(:,1),pairsOfLengthsNeighbors(:,2),'type','Pearson');
    end
    
    allpairsNforFrame = numel(pairsOfLengths); 
    allpairsN(end+1) = allpairsNforFrame;
    
    neighborNforFrame = numel(pairsOfLengthsNeighbors);
    neighborN(end+1) = neighborNforFrame;
    
    disp(['Selected ' num2str(neighborNforFrame) '/' num2str(allpairsNforFrame) ' as neighbors']);
    
    
    % Now select pairs that have certain distance
    
    
    %{
    std1 = std(pairsOfLengths(:,1));
    std2 = std(pairsOfLengths(:,2));
    mycorrManual = mean((pairsOfLengths(:,1).*pairsOfLengths(:,2)))/(std1*std2)
    %}

    disp(['Analyzing range ' num2str(min(myRange)) '-' num2str(max(myRange)) ', now at frame '  num2str(fr)]);
    
end

%% plot
figure (3); 
subplot(1,2,1); hold on;
plot([min(myRange),max(myRange)],[0,0],'k');
l1=plot(myRange,correlationsPerFrame,'k','Linewidth',2);
l2=plot(myRange,correlationsPerFrameNeighbors,'b','Linewidth',2);
ylim([-1,1]);
xlim([min(myRange),max(myRange)])
%legend(l1,'All pairs');

subplot(1,2,2); hold on;
l1=plot(myRange,allpairsN,'k--','Linewidth',3);
l2=plot(myRange,neighborN,'b.','Linewidth',3);
xlim([min(myRange),max(myRange)])
%legend(l1,'All pairs');

%% 
myFrame = 130;
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
