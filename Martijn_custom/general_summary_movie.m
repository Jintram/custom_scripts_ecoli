
% initiate
% myRange=[4:10]; % custom range
myRange = settings.frameRangeFull; 
croppingRecteangle = [p.cropLeftTop p.cropRightBottom];

% obtain paramter to plot
% length
[timesPerframe_nrs, paramPerframe_nrs, schnitzesPerFrame_nrs, plottingfXaxis] = MW_getparamfromschnitzcells(schnitzcells, 'time', 'length_fitNew');
% fluor signal
[fluorTimesPerframe_nrs, fluorParamPerframe_nrs, schnitzesPerFrame_nrs, fluorPlottingfXaxis] = MW_getparamfromschnitzcells(schnitzcells, 'time', 'G6_mean_all');

% generate some colors for plotting
totalNrcolors = max([schnitzesPerFrame_nrs{:}]);
preferredlinecolors = linspecer(totalNrcolors)
randomizationIndices = randperm(totalNrcolors);
preferredlinecolors(:,1) = preferredlinecolors(randomizationIndices,1)
preferredlinecolors(:,2) = preferredlinecolors(randomizationIndices,2)
preferredlinecolors(:,3) = preferredlinecolors(randomizationIndices,3)


% set up plots
h = figure(2); clf; hold on;
% plot for length 
subplot(2,2,2); 
myXlim = [min([timesPerframe_nrs{:}]) max([timesPerframe_nrs{:}])];
myYlim = [min([paramPerframe_nrs{:}]) max([paramPerframe_nrs{:}])];
xlim(myXlim);
ylim(myYlim);
% plot for fluor 
subplot(2,2,4); 
fluormyXlim = [min([fluorTimesPerframe_nrs{:}]) max([fluorTimesPerframe_nrs{:}])];
fluormyYlim = [min([fluorParamPerframe_nrs{:}]) max([fluorParamPerframe_nrs{:}])];
xlim(fluormyXlim);
ylim(fluormyYlim);

% loop over the frames
for fr = myRange
    
  % load the seg file (contains segmented and phase img)
  name= [p.segmentationDir,p.movieName,'seg',str3(fr)];
  load(name);
  
  % Create colony picture using PN_imshowlabel and plot
  subplottight(2,2,1);
  p.showPerim = 1;
  outimage = PN_imshowlabel(p, Lc, croppingRecteangle,[],[],'phaseImage',phsub); 
  imshow(outimage,[])
  %text(20,20,['frame ' num2str(fr)],'Color','k','FontWeight','bold','BackgroundColor','white');

  % Colony image, segmented
  subplottight(2,2,3);
  p.showPerim = 0;
  outimage = PN_imshowlabel(p, Lc, croppingRecteangle,[],[]);
  imshow(outimage,[])
  %text(20,20,['frame ' num2str(fr)],'Color','k','FontWeight','bold','BackgroundColor','white');  
  
  % Plot growth of schnitzcells
  subplot(2,2,2); hold on;
  % plot each schnitz separate color
  for i=1:numel(schnitzesPerFrame_nrs{fr})
    plot(timesPerframe_nrs{fr}(i), paramPerframe_nrs{fr}(i),'.',...
           'MarkerSize',10,'LineWidth',3,'MarkerEdgeColor',preferredlinecolors(schnitzesPerFrame_nrs{fr}(i),:));  
  end
  % Cosmetics
  MW_makeplotlookbetter(20);
  xlabel('Time [mins]')
  ylabel('Bacterial length [um]')
  title(['Frame nr ' num2str(fr)])
  
  % Plot fluor of schnitzcells
  subplot(2,2,4); hold on;
  % plot each schnitz separate color
  for i=1:numel(schnitzesPerFrame_nrs{fr})
        plot(fluorTimesPerframe_nrs{fr}(i), fluorParamPerframe_nrs{fr}(i),'.',...
            'MarkerSize',10,'LineWidth',3,'MarkerEdgeColor',preferredlinecolors(schnitzesPerFrame_nrs{fr}(i),:));   
  end
  % Cosmetics
  MW_makeplotlookbetter(20);
  xlabel('Time [mins]')
  ylabel('Fluor signal [a.u.]')
  
  % save figure
  saveas(h, ['D:\Local_Playground\mymovietestv2\movietest_' sprintf('%04d',fr) '.jpg']);    
  
end



%{
figure(1), imshow(Lc,[])
figure(2), imshow(phsub,[])

figure(2), imshow(outim,[])
%}
