
p = DJK_initschnitz('pos8crop','2014-05-01','e.coli.AMOLF','rootDir','F:\A_Tans1_step4a_partially_analyzed_analysis\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
myRange=[39:130]

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
% parents.
% First pick some colours for the young guns:
nrYoungGuns = numel(youngSchnitzes);
YoungGunColors = distinguishable_colors(nrYoungGuns);
% Now label each schnitz for its relatedness to the final offspring
IncestColor = cell(size(schnitzcells));
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
        IncestColor{i} = mean(allYoungColors);
    else
        IncestColor{i} = allYoungColors;
    end;
end

% figure
h=figure(1), set(h, 'Position', [0,0,800,500]);

% loop over the frames
for fr = myRange
    
  % Plot the image
  % ===
    
  % load the seg file (contains segmented and phase img)
  name= [p.segmentationDir,p.movieName,'seg',str3(fr)];
  load(name);
  
  % 
  p.showPerim = 1;
  outim = PN_imshowlabel(p, Lc,0,0,0,'phaseImage',phsub); % note i'm feeding the same img as "previous" img
  
  % crop img
  %outim = imcrop(outim, [395   225   650   650]);
  
  h=figure(1), subplottight(1,2,1)
  imshow(outim,[])
  text(20,20,['frame ' sprintf('%05d', fr)],'Color','k','FontWeight','bold','BackgroundColor','white');  
  
  % Now also plot length
  % ===
  
  [theSchnitzes,cellTimePoints,cellLengths,timemin,timemax,lengthmin,lengthmax] =  [theSchnitzes,cellTimePoints,cellLengths,timemin,timemax,lengthmin,lengthmax] =  givemeschnitzinfoforframe(p, schnitzcells, fr);
  
  subplot(1,2,2), hold on 
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

%{
figure(1), imshow(Lc,[])
figure(2), imshow(phsub,[])

figure(2), imshow(outim,[])
%}
