function h = timeevolutionplot(time,y1,y2,highlightpoints,highlightPointCategories,figNr)

categorieColors= linspecer(max(highlightPointCategories)+1);

h=figure(figNr); clf; hold on;
title('Time evolution');

% plot
lastIdx = numel(time)-1;
timeColors = makeColorMap([0 0 1], [0 0 0], lastIdx);
timeColors = colormap(parula(lastIdx));
%arrowPlot(CRPvalues,growthvalues,'number',10,'color',[0 0 0]); hold on;
switchcount=0; switchH = [];
for timeIdx = 1:lastIdx    
    plot([y1(timeIdx) y1(timeIdx+1)],[y2(timeIdx) y2(timeIdx+1)],'x-','Color',timeColors(timeIdx,:),'LineWidth',2)
    
    if ismember(timeIdx, highlightpoints)
        switchcount=switchcount+1;
        switchColor = categorieColors(highlightPointCategories(switchcount)+1,:);
        switchH(end+1)=plot(y1(timeIdx),y2(timeIdx),'s','MarkerSize',7,'LineWidth',3,'Color',switchColor); %'MarkerFaceColor','k',
    end
    
    arrowu=y1(timeIdx+1) -y1(timeIdx);
    arrowv=y2(timeIdx+1) -y2(timeIdx);
    %quiver(CRPvalues(timeIdx),growthvalues(timeIdx),arrowu,arrowv);
    %arrow([CRPvalues(timeIdx) growthvalues(timeIdx)],[CRPvalues(timeIdx+1) growthvalues(timeIdx+1)],'Length',15,'Width',2);
end

legend(switchH(1:2),{'2100uM cAMP','43uM'});

xlabel('cAMP.CRP reporter');
ylabel('Growth rate');

colormap(timeColors);
hC=colorbar();

ylabel(hC, 'Time (hrs)')

inputSettings.rangeIn = [0,1]; % original range of axis
inputSettings.desiredSpacing = 2; %  desired spacing of axis in new metric
inputSettings.rangeOut = [0,max(time)]; % the desired target range
[tickLocationsOldMetric, correspondingLabels] = labelremapping(inputSettings);

set(hC, 'YTickLabel',correspondingLabels, 'XTick',tickLocationsOldMetric)

MW_makeplotlookbetter(20);

