

%% Analyze data
meanWidth={};
allWidth=[];
flattenedallLengthsOfBacteriaInMicronsMW=[];
for dataSetIdx = 1:5

    % Load dataset
    switch dataSetIdx
        case 1
            % Dataset 1
            %load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-09\pos3crop\analysis\straightenedCells\2013-12-09pos3crop_straightFluorData.mat');
            load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-09\pos3crop\data\pos3crop-skeletonData.mat');
        case 2
            % Dataset 2
            load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos4crop\analysis\straightenedCells\2013-09-24pos4crop_straightFluorData.mat');
            load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos4crop\data\pos4crop-skeletonData.mat');
        case 3
            % Dataset 3
            load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos4crop\analysis\straightenedCells\2013-12-16pos4crop_straightFluorData.mat');
            load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos4crop\data\pos4crop-skeletonData.mat');
        case 4
            % Dataset 4
            load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos5crop\analysis\straightenedCells\2013-09-24pos5crop_straightFluorData.mat');
            load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos5crop\data\pos5crop-skeletonData.mat');
        case 5
            % Dataset 5
            load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos5crop\analysis\straightenedCells\2013-12-16pos5crop_straightFluorData.mat');
            load('G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos5crop\data\pos5crop-skeletonData.mat');
    end
    
    %% analyze it
    for frIdx = 1:numel(allbacterialWidthMicrons)
        if ~isempty(allbacterialWidthMicrons{frIdx})

            for cellIdx = 1:numel(allbacterialWidthMicrons{frIdx})

                %% 
                bacNPixels = numel(allbacterialWidthMicrons{frIdx}{cellIdx});
                selectedIdxs = [round(bacNPixels.*.2):round(bacNPixels*.4),...
                                round(bacNPixels.*.6):round(bacNPixels*.8)];
                meanWidth{frIdx}{cellIdx} = mean(allbacterialWidthMicrons{frIdx}{cellIdx}(selectedIdxs));
                allWidth(end+1)=mean(allbacterialWidthMicrons{frIdx}{cellIdx}(selectedIdxs));

                %% also save calculated lengths
                flattenedallLengthsOfBacteriaInMicronsMW(end+1) = allLengthsOfBacteriaInMicrons{frIdx}(cellIdx);

            end

        end
    end
    
    disp(['Analysis #' num2str(dataSetIdx) 'finished']);
    
end
%%
figure
[counts,edges]=histcounts(allWidth,30);
binCenters = [edges(2:end)+edges(1:end-1)]./2;
plot(binCenters,counts);

%%
figure
[counts,edges]=histcounts(flattenedallLengthsOfBacteriaInMicronsMW,30);
binCenters = [edges(2:end)+edges(1:end-1)]./2;
plot(binCenters,counts);


%% get some stats
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins]=...
    binnedaveraging({flattenedallLengthsOfBacteriaInMicronsMW},{allWidth},[0:3:40]);

%%
hBacWidths=figure; clf; hold on;
scatter(flattenedallLengthsOfBacteriaInMicronsMW,allWidth,...
     7^2,[.5 .5 .5],'filled','MarkerFaceAlpha',.5);
scatter(binCenters,meanValuesForBins,...
    15^2,'k','filled');
xlabel('Bacterial Length (µm)');
ylabel('Bacterial Width (µm)');
MW_makeplotlookbetter(15);
xlim([0,40]); ylim([0,1]);
%title(['N=' num2str(numel(allWidth))]);

% bin info
%myColor=[255 187 84]./255;
myColor=[0 0 0];
for i=1:numel(counts)
    if counts(i)>0
        text(binCenters(i),.01,...
            ['N=' num2str(counts(i))],...
            'Color',myColor,'FontSize',12,'FontWeight','bold','Rotation',90)
    end
end
plot([0:3:50],zeros(numel([0:3:50]),1)','^','MarkerFaceColor',myColor,'MarkerEdgeColor',myColor)

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

%%

%binnedpdfs %binnedaveraging
hWidthViolin=figure;
violin({binnedValues{1:end-2}},'x',binCenters)
xlim([0,30]); ylim([0,1]);

% bin info
%myColor=[255 187 84]./255;
myColor=[0 0 0];
for i=1:7
    if counts(i)>0
        text(binCenters(i),.01,...
            ['N=' num2str(counts(i))],...
            'Color',myColor,'FontSize',12,'FontWeight','bold','Rotation',90)
    end
end
plot([0:3:50],zeros(numel([0:3:50]),1)','^','MarkerFaceColor',myColor,'MarkerEdgeColor',myColor)

%cosmetics
xlabel('Bacterial Length (µm)');
ylabel('Bacterial Width (µm)');
MW_makeplotlookbetter(15);
xlim([0,21]); ylim([0,1]);


%%
figure
boxplot(binnedValues{1})


