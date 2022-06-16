
%%
allSizes=[];
for idx=1:numel(datasetsPaths)
    
    schnitzcells = loadandrename(datasetsPaths{idx});
    cellSizes=arrayfun(@(i) schnitzcells(i).length_fitNew(1), 1:numel(schnitzcells));
    allSizes=[allSizes cellSizes];
    disp('.');
    
end
  
%%
figure; clf; hold on;
[n,e]=histcounts(allSizes,100);
c=e(2:end)-(e(2)-e(1))/2;
plot(c,n,'-r','LineWidth',2)
plot(1.5.*1.8.*[1,3,5,8],[0,0,0,0],'^','MarkerFaceColor','k','MarkerSize',7);
ylim([0,200]);
xlabel('Cell birth size');
ylabel('Count');
MW_makeplotlookbetter(20);