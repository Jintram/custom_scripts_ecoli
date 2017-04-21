
% Set up
clear combinedDataAdder;
YLIMbirthlife=[0,150];

%% First collect two other conditions

% SulA
WHATDATA = 'sulA';
NOSAVEPLEASE=1;
LENGTHFIELD = 'length_skeleton';

RUNSECTIONSFILADIV = 'loadData'; 
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'adderLengthLifetime';
script20160429_filamentRecoveryDivisionRatioss

%% combine data
combinedDataAdder.('sulA').alldatax = alldatax;
combinedDataAdder.('sulA').alldatay = alldatay;
combinedDataAdder.('sulA').meanValuesForBins = meanValuesForBins;
combinedDataAdder.('sulA').binCenters = binCenters;
combinedDataAdder.('sulA').stdValuesForBins = stdValuesForBins;

%% Temperature
WHATDATA = 'temperature';
NOSAVEPLEASE=1;
LENGTHFIELD = 'length_skeleton';

RUNSECTIONSFILADIV = 'loadData'; 
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'adderLengthLifetime';
script20160429_filamentRecoveryDivisionRatioss

%% combine data
combinedDataAdder.('temperature').alldatax = alldatax;
combinedDataAdder.('temperature').alldatay = alldatay;
combinedDataAdder.('temperature').meanValuesForBins = meanValuesForBins;
combinedDataAdder.('temperature').binCenters = binCenters;
combinedDataAdder.('temperature').stdValuesForBins = stdValuesForBins;


%%
WHATDATA = 'tetracycline';
NOSAVEPLEASE=1;
LENGTHFIELD = 'length_fitNew';

YLIMbirthlife=[0,200];

RUNSECTIONSFILADIV = 'loadData'; % (same as a but reload when executed per section)
script20160429_filamentRecoveryDivisionRatioss

RUNSECTIONSFILADIV = 'adderLengthLifetime';
script20160429_filamentRecoveryDivisionRatioss

%% combine data
combinedDataAdder.('tetracycline').alldatax = alldatax;
combinedDataAdder.('tetracycline').alldatay = alldatay;
combinedDataAdder.('tetracycline').meanValuesForBins = meanValuesForBins;
combinedDataAdder.('tetracycline').binCenters = binCenters;
combinedDataAdder.('tetracycline').stdValuesForBins = stdValuesForBins;

%%
DATASETNAMES={'tetracycline','temperature','sulA'};
LEGENDNAMES = {'Tetracycline','Temperature','SulA'};

somecolors=linspecer(3);
mycolors = [0, 0, 0; somecolors(2:3,:)];

hAdderFigs=[];
for dataIdx = 1:3
    hAdderFigs(dataIdx) = figure(dataIdx); clf; hold on;
    
    x=combinedDataAdder.(DATASETNAMES{dataIdx}).alldatax;
    y=combinedDataAdder.(DATASETNAMES{dataIdx}).alldatay;
    scatter(x,y,7^2,...mycolors(dataIdx,:),...
        [.5 .5 .5],'filled','MarkerFaceAlpha',.3);
%end
%for dataIdx = 1:3
    
    y=combinedDataAdder.(DATASETNAMES{dataIdx}).meanValuesForBins;
    x=combinedDataAdder.(DATASETNAMES{dataIdx}).binCenters;
    err=combinedDataAdder.(DATASETNAMES{dataIdx}).stdValuesForBins;

    lErr=errorbar(x,y,err,'o-','LineWidth',2,'Color',mycolors(dataIdx,:))
   
%end
%for dataIdx = 1:3
    
    y=combinedDataAdder.(DATASETNAMES{dataIdx}).meanValuesForBins;
    x=combinedDataAdder.(DATASETNAMES{dataIdx}).binCenters;

    plot(x,y,'o-','LineWidth',2,'Color',mycolors(dataIdx,:),'MarkerFaceColor',mycolors(dataIdx,:))
   
    %
    xlim([0,40]);
    if dataIdx==1, xlim([0,20]); end
    ylim([0,6]);

    % 
    MW_makeplotlookbetter(20);
    xlabel('Birth size');
    ylabel('Added length (\mum)');
    
    % 
    legend(lErr,LEGENDNAMES{dataIdx});
    
end


    
    
    
    
    