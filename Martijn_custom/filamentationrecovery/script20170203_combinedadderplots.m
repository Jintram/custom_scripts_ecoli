
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
combinedDataAdder.('sulA').stdErrValuesForBins = stdErrValuesForBins;

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
combinedDataAdder.('temperature').stdErrValuesForBins = stdErrValuesForBins;


%% Then tetracycline
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
combinedDataAdder.('tetracycline').stdErrValuesForBins = stdErrValuesForBins;

%% Also do the deltamin one

% deltaMinTET
WHATDATA = 'deltaMinTET';
NOSAVEPLEASE=1;
LENGTHFIELD = 'length_skeleton';

RUNSECTIONSFILADIV = 'loadData'; 
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'adderLengthLifetime';
script20160429_filamentRecoveryDivisionRatioss

% combine data
combinedDataAdder.(WHATDATA).alldatax = alldatax;
combinedDataAdder.(WHATDATA).alldatay = alldatay;
combinedDataAdder.(WHATDATA).meanValuesForBins = meanValuesForBins;
combinedDataAdder.(WHATDATA).binCenters = binCenters;
combinedDataAdder.(WHATDATA).stdValuesForBins = stdValuesForBins;
combinedDataAdder.(WHATDATA).stdErrValuesForBins = stdErrValuesForBins;

%%
DATASETNAMES={'tetracycline','temperature','sulA','deltaMinTET'};
LEGENDNAMES = {'Tetracycline','Temperature','SulA','deltaMinTET'};

somecolors=linspecer(5);
mycolors = [0, 0, 0; somecolors(2:4,:)];

hAdderFigs=[];
for dataIdx = 1:4
    hAdderFigs(dataIdx) = figure(dataIdx); clf; hold on;
    
    x=combinedDataAdder.(DATASETNAMES{dataIdx}).alldatax;
    y=combinedDataAdder.(DATASETNAMES{dataIdx}).alldatay;
    scatter(x,y,7^2,...mycolors(dataIdx,:),...
        [.5 .5 .5],'filled','MarkerFaceAlpha',.3);
%end
%for dataIdx = 1:3
    
    y=combinedDataAdder.(DATASETNAMES{dataIdx}).meanValuesForBins;
    x=combinedDataAdder.(DATASETNAMES{dataIdx}).binCenters;
    err=combinedDataAdder.(DATASETNAMES{dataIdx}).stdErrValuesForBins;

    lErr=errorbar(x,y,err,...                
                '-o','Color',mycolors(dataIdx,:),'MarkerSize',5,'LineWidth',2)
                %'o-','LineWidth',2,'Color',mycolors(dataIdx,:))
   
%end
%for dataIdx = 1:3
    
    y=combinedDataAdder.(DATASETNAMES{dataIdx}).meanValuesForBins;
    x=combinedDataAdder.(DATASETNAMES{dataIdx}).binCenters;

    %plot(x,y,...
    %    '-o','Color',mycolors(dataIdx,:),'MarkerSize',5,'LineWidth',2);
        %'o-','LineWidth',2,'Color',mycolors(dataIdx,:),'MarkerFaceColor',mycolors(dataIdx,:))
   
    %
    xlim([0,40]);
    if dataIdx==1, xlim([0,20]); end
    ylim([0,6]);
    if dataIdx==4, ylim([0,13]); end

    % 
    MW_makeplotlookbetter(20);
    xlabel('Birth size (µm)');
    ylabel('Added length (µm)');
    
    % Legend
    % ==
    %legend(lErr,LEGENDNAMES{dataIdx});
    %legend boxoff
    
end


    
    
    
    
    