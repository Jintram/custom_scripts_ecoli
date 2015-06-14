
p = DJK_initschnitz('pos8crop','2014-05-01','e.coli.AMOLF','rootDir','F:\A_Tans1_step4a_partially_analyzed_analysis\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
myRange=[133:150]

% For plotting lengths
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

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
  outim = PN_imshowlabel(p, Lc,[],[],[],'phaseImage',phsub); % note i'm feeding the same img as "previous" img
  
  % crop img
  %outim = imcrop(outim, [395   225   650   650]);
  
  h=figure(1), subplottight(1,2,1)
  imshow(outim,[])
  text(20,20,['frame ' sprintf('%05d', fr)],'Color','k','FontWeight','bold','BackgroundColor','white');  
  
  % Now also plot length
  % ===
  
  [theSchnitzes,cellTimePoints,cellLengths,timemin,timemax,lengthmin,lengthmax] =  givemeschnitzinfoforframe(p, schnitzcells, fr);
  
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
