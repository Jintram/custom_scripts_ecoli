
%% Settings
BINEDGES = [1:1:40]; %BINEDGES = [0:3:40];
XCUTOFFVIOLIN=21;
XCUTOFFSCATTER=30;

%% In case you want to use the copied datasets. 
%{
WHATDATA='tetracycline';
RUNSECTIONSFILADIV='none';
NOSAVEPLEASE=1;
script20160429_filamentRecoveryDivisionRatioss
%}

%% updatedDatasets

% These are the original datasets, which have been updated with analyses of
% skeleton length.
originaldatasetsPaths={...
    'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-09\pos3crop\data\pos3crop-Schnitz.mat', ...
    'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos4crop\data\pos4crop-Schnitz.mat',...
    'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos4crop\data\pos4crop-Schnitz.mat',...
    'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos5crop\data\pos5crop-Schnitz.mat',...
    'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos5crop\data\pos5crop-Schnitz.mat',...
};
            
%%


allDatasetsWidths = [];
allDatasetsLengths = [];
allDatasetsColors =[];
LENGTHFIELD = 'length_fitNew';
LENGTHFIELD = 'length_skeleton';
for dataIdx = 1:4%numel(originaldatasetsPaths)
            
    %% 
    % load data
    %schnitzcells = loadandrename(datasetsPaths{dataIdx});
    load(originaldatasetsPaths{dataIdx});
    
    % simply brute force this
    allWidths=[]; allLengths=[];
    myColors = linspecer(numel(schnitzcells));
    colorsMixed = myColors(randperm(size(myColors,1)),1);
    colorsMixed = myColors(randperm(size(myColors,1)),2);
    colorsMixed = myColors(randperm(size(myColors,1)),3);
    allColors=[];
    
    for schnitzIdx = 1:numel(schnitzcells)

        %%
        lenghts = schnitzcells(schnitzIdx).(LENGTHFIELD);
        areas   = schnitzcells(schnitzIdx).area;
        widths  = areas./lenghts;

        
        allWidths  = [allWidths, widths];
        allLengths = [allLengths, lenghts];
        allColors  = [allColors; repmat(colorsMixed(schnitzIdx,:), [numel(lenghts),1])];

    end

    allDatasetsWidths = [allDatasetsWidths, allWidths];
    allDatasetsLengths = [allDatasetsLengths, allLengths];
    allDatasetsColors = [allDatasetsColors; allColors];
    
    disp(['Analysis for dataset # ' num2str(dataIdx) ' done.']);
end


%% simple figure
figure; clf; hold on;
scatter(allDatasetsLengths,allDatasetsWidths,5^2,allDatasetsColors); 

[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins] = ...
    binnedaveraging({allDatasetsLengths},{allDatasetsWidths},[0:3:35]);

scatter(binCenters,medianValuesForBins,20^2,'k','FaceColor','k')


%%
figure
[counts,edges]=histcounts(allDatasetsWidths,30);
binCenters = [edges(2:end)+edges(1:end-1)]./2;
plot(binCenters,counts);

%%
figure
[counts,edges]=histcounts(allDatasetsLengths,30);
binCenters = [edges(2:end)+edges(1:end-1)]./2;
plot(binCenters,counts);


%% get some stats
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins]=...
    binnedaveraging({allDatasetsLengths},{allDatasetsWidths},BINEDGES);

%%
hBacWidths=figure(4); clf; hold on;
scatter(allDatasetsLengths,allDatasetsWidths,...
     3^2,[.5 .5 .5],'filled','MarkerFaceAlpha',.5);
scatter(binCenters,meanValuesForBins,...
    7^2,'k','filled');
xlabel('Bacterial Length (µm)');
ylabel('Bacterial Width (µm)');
MW_makeplotlookbetter(15);
xlim([0,XCUTOFFSCATTER]); ylim([0,1]);
%title(['N=' num2str(numel(allDatasetsWidths))]);

% bin info
%myColor=[255 187 84]./255;
myColor=[0 0 0];
for i=1:numel(BINEDGES(BINEDGES<XCUTOFFSCATTER))
    if counts(i)>0
        text(binCenters(i),.01,...
            ['N=' num2str(counts(i))],...
            'Color',myColor,'FontSize',8,'FontWeight','bold','Rotation',90)
    end
end
plot(BINEDGES,zeros(numel(BINEDGES),1)','^','MarkerFaceColor',myColor,'MarkerEdgeColor',myColor,'MarkerSize',4)

figure(4);hold on;
%linspecer(3)
%W=.72;
W=.74;
CAPSIZE=W^2/2;
CAPSIZE=W^2-(W/2)^2*pi;
plot([[1:30].*W-CAPSIZE]./[1:30],'LineWidth',2,'Color',[0.9153    0.2816    0.2878]);
ylim([0,1]);

%%
disp('Two sample t-test');
disp('===');
for otherbin = 1:10
    COMPAREBINS=[1,otherbin];
    [h,p] = ttest2(binnedValues{COMPAREBINS(1)},binnedValues{COMPAREBINS(2)});
    % returns a test decision for the null hypothesis that the data in vectors
    % x and y comes from independent random samples from normal distributions
    % with equal means and equal but unknown variances
    % 0 = does not reject [same normal distr.]
    % 1 = rejects null-hypothesis 

    if h
        disp(['H_0 X rejected. Bins ' num2str(COMPAREBINS(1)) ' and ' num2str(COMPAREBINS(2)) ' according to two sample t-test are not from same distribution (p=' num2str(p) ').']);
    else
        disp(['H_0 A accepted. Bins ' num2str(COMPAREBINS(1)) ' and ' num2str(COMPAREBINS(2)) ' according to two sample t-test are from same distribution (p=' num2str(p) ').']);
    end
end

%% Further stats
 percentageWidths = meanValuesForBins./meanValuesForBins(1)
 meanFirstBin = mean(meanValuesForBins(1))

%% Violin plots

%binnedpdfs %binnedaveraging
hWidthViolin=figure;
violin({binnedValues{1:end-2}},'x',binCenters)
xlim([0,30]); ylim([0,1]);

% bin info
%myColor=[255 187 84]./255;
myColor=[0 0 0];
for i=1:numel(BINEDGES(BINEDGES<XCUTOFFVIOLIN))
    if counts(i)>0
        text(binCenters(i),.01,...
            ['N=' num2str(counts(i))],...
            'Color',myColor,'FontSize',12,'FontWeight','bold','Rotation',90)
    end
end
plot(BINEDGES(BINEDGES<XCUTOFFVIOLIN),zeros(numel(BINEDGES(BINEDGES<XCUTOFFVIOLIN)),1)','^','MarkerFaceColor',myColor,'MarkerEdgeColor',myColor)

%cosmetics
xlabel('Bacterial Length (µm)');
ylabel('Bacterial Width (µm)');
MW_makeplotlookbetter(15);
xlim([0,XCUTOFFVIOLIN]); ylim([0,1]);


%%
figure
boxplot(binnedValues{1})


