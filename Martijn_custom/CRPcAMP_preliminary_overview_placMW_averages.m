
avgsY=[]; avgsC=[];
stdsY=[]; stdsC=[];
some_colors;

myGrouping = [ 1,2,3,4 ...
   ];

newDataSets = { ...
    'F:\A_Tans1_step1_incoming_not_backed_up\2015-04-04\553ace50mscrop\data\553ace50mscrop-Schnitz.mat', ...    
    'F:\A_Tans1_step1_incoming_not_backed_up\2015-04-04\553lac50mscrop\data\553lac50mscrop-Schnitz.mat', ...
    'F:\A_Tans1_step1_incoming_not_backed_up\2015-04-04\553mal50mscrop\data\553mal50mscrop-Schnitz.mat' ...
    'F:\A_Tans1_step1_incoming_not_backed_up\2015-04-04\553RDM50mscrop\data\553RDM50mscrop-Schnitz.mat' ...
    }
    
% TODO annotate better
% noreen's mu from ace, lac, mal, RDM
noreenmu = [0.195012294777113, 0.5873, 0.6209, 1.7809];


%% final figure w. averages Y
% -----------
figure(7), clf; hold on;

legendLines = [];
for i = [1:4]
    load(newDataSets{i});
    myConcentrations = [schnitzcells.av_Y6_mean];
    mus = ones(numel(myConcentrations))*noreenmu(i);
    plot(mus,myConcentrations,'x','MarkerSize',10,'MarkerEdgeColor','k')    
    
    theMedian = median(myConcentrations);
    theStd = std(myConcentrations);
    %errorbar(noreenmu(i),theMedian,theStd,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15);
    lineH = plot(noreenmu(i),theMedian,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15);
    legendLines = [legendLines lineH];
end
legend(legendLines, {'ace', 'lac', 'mal', 'RDM'})

ylim([0 max(myConcentrations)*2 ]);

xlabel('Growth rate (dbl/hr) [exp NW]');
ylabel('Concentration [my exp]');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);



title('YFP data')

% Plot w. C data
% ------------------
figure(8), clf; hold on;

legendLines = [];
for i = [1:4]
    load(newDataSets{i});
    myConcentrations = [schnitzcells.av_C6_mean];
    mus = ones(numel(myConcentrations))*noreenmu(i);
    plot(mus,myConcentrations,'x','MarkerSize',10,'MarkerEdgeColor','k')    
    
    theMedian = median(myConcentrations);
    theStd = std(myConcentrations);
    %errorbar(noreenmu(i),theMedian,theStd,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15);
    lineH = plot(noreenmu(i),theMedian,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15);
    legendLines = [legendLines lineH];
end
legend(legendLines, {'ace', 'lac', 'mal', 'RDM'})

ylim([0 max(myConcentrations)*2 ]);

xlabel('Growth rate (dbl/hr) [exp NW]');
ylabel('Concentration [my exp]');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

title('CFP data')




