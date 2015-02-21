
% Settings
% ===
ANCESTORDILUTION = .75; % How much color mixing dies out from subsequent ancestry
myRange=[39:130]

% Load schnitz stuff
p = DJK_initschnitz('pos8crop','2014-05-01','e.coli.AMOLF','rootDir','F:\A_Tans1_step4a_partially_analyzed_analysis\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');

% For plotting lengths
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

% THIS CODE DOES EXACTLY THE OPOSITE!
%{ 
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
%}
aColorMap = distinguishable_colors(numel(schnitzcells), [0,0,0]);


% figure
h=figure(1), set(h, 'Position', [0,0,800,500]);

% loop over the frames
h=figure(1); 
for fr = myRange
    
  % Plot the image
  % ===
    
  % load the seg file (contains segmented and phase img)
  name= [p.segmentationDir,p.movieName,'seg',str3(fr)];
  load(name);        
   
  % Get some info on this frame
  % ===
  [theSchnitzes,cellTimePoints,cellLengths,timemin,timemax,lengthmin,lengthmax,mapIndexes] = ...
      givemeschnitzinfoforframe(p, schnitzcells, fr);
  
  % Now we edit color map such that sisters look like their mothers, and a
  % bit less like their grandmothers, etc.
  ancestryMap = ones(size(schnitzcells,2)+1,3);
  for i = theSchnitzes
      colorsAncestry = []; currentAncestor = i; done = 0;
      relatedness = 1;
      while ~done
          colorsAncestry(end+1,:) = aColorMap(currentAncestor,:).*relatedness;
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
  for i = 1:numel(mapIndexes)
      myCustomColorMap(mapIndexes(i)+1,:) = ancestryMap(theSchnitzes(i),:);
  end    
    
  % make image of cells
  %p.showPerim = 1;
  outim = PN_imshowlabel(p, Lc,0,0,0,'customColors',myCustomColorMap);%,'phaseImage',); % note i'm feeding the same img as "previous" img
  % crop img
  %outim = imcrop(outim, [395   225   650   650]);
  
  % Plot image of cells and length plot
  % ===
  
  %  image of cells
  subplottight(1,2,1);
  imshow(outim,[]);
  text(20,20,['frame ' sprintf('%05d', fr)],'Color','k','FontWeight','bold','BackgroundColor','white');  
  
  % their lengths in a plot
  subplot(1,2,2), hold on ;
  for j = 1:numel(theSchnitzes)
      
    %plot(ones(numel(cellLengths),1).*fr, cellLengths, 'o','Color',ind2rgb(j,cmap))
    plot(fr, cellLengths(j), '.','Color',ind2rgb(j,colormap(lines)))
  end
  text(20,20,['frame ' num2str(fr)],'Color','k','FontWeight','bold','BackgroundColor','white');
  xlim([min(myRange),max(myRange)]), ylim([lengthmin,lengthmax]); % redundant  
  
  title('Cell lengths','FontSize',15)
  xlabel('Frame number','FontSize',15)
  ylabel('Cell length ({\mu}m)','FontSize',15)
  
  
  
  saveas(h, ['D:\Local_Playground\mymovietest\movietest_lengths_' sprintf('%05d',fr) '.jpg']);
  
end

msgbox('Done');

%{
figure(1), imshow(Lc,[])
figure(2), imshow(phsub,[])

figure(2), imshow(outim,[])
%}
