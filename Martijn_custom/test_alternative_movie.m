
p = DJK_initschnitz('pos1crop','2014-05-01','e.coli.AMOLF','rootDir','F:\A_Tans1_step4a_partially_analyzed_analysis\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
myRange=[39:130]

% loop over the frames
for fr = myRange
    
  % load the seg file (contains segmented and phase img)
  name= [p.segmentationDir,p.movieName,'seg',str3(fr)];
  load(name);
  
  % 
  p.showPerim = 1;
  outim = PN_imshowlabel(p, Lc,rect,Lc,rect,'phaseImage',phsub); % note i'm feeding the same img as "previous" img
  
  % crop img
  outim = imcrop(outim, [395   225   650   650]);
  
  h=figure(2), imshow(outim,[])
  text(20,20,['frame ' num2str(fr)],'Color','k','FontWeight','bold','BackgroundColor','white');
  saveas(h, ['D:\Local_Playground\mymovietest\movietest_' num2str(fr) '.jpg']);
  
end

%{
figure(1), imshow(Lc,[])
figure(2), imshow(phsub,[])

figure(2), imshow(outim,[])
%}
