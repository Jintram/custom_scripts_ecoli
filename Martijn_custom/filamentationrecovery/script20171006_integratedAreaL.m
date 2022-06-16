

%% define output paths

OUTPUTDIR = 'D:\Local_Data\Dropbox\Dropbox\Filamentation recovery\Submission_Current_Biology_round2\RebuttalFigs\matlabExport\'


%% Load some info about the data (e.g. paths)

NOSAVEPLEASE=1;
RUNSECTIONSFILADIV='none';
script20160429_filamentRecoveryDivisionRatioss;


%% Calculating integrated area under bacterial size
allAreas          = [];
%allAreasPrime     = [];
allMotherLengths  = [];
allAverageMus     = [];
for dataSetIndex = 1:numel(datasetsPaths)

    % Load schnitz
    if ~exist('simulatedschnitzcells','var')
        if ~SPECIALCASE
            schnitzcells = loadandrename(datasetsPaths{dataSetIndex});
        else
            schnitzcells=S_all_shifted{dataSetIndex};
        end
    end
        
    % Let's simply brute force this
    Areas=NaN(1, numel(schnitzcells));
    %AreasPrime=NaN(1, numel(schnitzcells));
    motherLengths=NaN(1, numel(schnitzcells));
    averageMu=NaN(1, numel(schnitzcells));
    for schnitzIdx = 1:numel(schnitzcells)

        %%     
        % get the data
        lRaw = schnitzcells(schnitzIdx).length_fitNew;
        %lRaw = schnitzcells(schnitzIdx).areaPixels;
        tRaw = schnitzcells(schnitzIdx).time;

        % perform the integration
        dt = (tRaw(2:end)-tRaw(1:end-1))/60; % divide by 60 to get hrs
        y  = (lRaw(2:end)+lRaw(1:end-1))/2;
        A  = sum(y.*dt);

        % perform integration of log(length)
        %Aprime  = sum(log(y).*dt);

        % collect data
        Areas(schnitzIdx)         = A;
        %AreasPrime(schnitzIdx)    = Aprime;
        motherLengths(schnitzIdx) = lRaw(1);
        
        averageMu(schnitzIdx)     = schnitzcells(schnitzIdx).av_mu_fitNew;

    end
    
    allAreas          = [allAreas Areas];
    %allAreasPrime     = [allAreasPrime AreasPrime];
    allMotherLengths  = [allMotherLengths motherLengths];
    
    allAverageMus     = [allAverageMus averageMu];
    
    disp(['Loop ' num2str(dataSetIndex) ' done']);
end

%%

%selectedDataIdx = allAreas<1000;
selectedDataIdx = 1:numel(allAreas);
disp('Note that we are throwing away some outliers here!');

hA=figure; clf; hold on;

[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins]=binnedaveraging({allMotherLengths(selectedDataIdx)},{allAreas(selectedDataIdx)},[0:2:40]);
L = binCenters;
A = medianValuesForBins;

scatter(allMotherLengths(selectedDataIdx),allAreas(selectedDataIdx),...
    5^2,[.7 .7 .7],'filled','MarkerFaceAlpha',.5);
%scatter(binCenters,meanValuesForBins,15^2,'r','filled')
scatter(L,A,7^2,'k','filled')

xlim([0,40]); 
ylim([0,20]);

%ylabel('Integrated L(t)dt (min*µm)');
ylabel(['Area under L(t) curve ' 10 ' (min*µm)']);
xlabel('Birth length (µm)');

MW_makeplotlookbetter(20);

%% Now also make a plot for the growth rates of the cells


hMu=figure; clf; hold on;

[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins]=binnedaveraging({allMotherLengths(selectedDataIdx)},{allAverageMus(selectedDataIdx)},[0:2:40]);
Lprime = binCenters;
growthRatesMean = medianValuesForBins;


scatter(allMotherLengths(selectedDataIdx),allAverageMus(selectedDataIdx),...
    5^2,[.7 .7 .7],'filled','MarkerFaceAlpha',.5);
%scatter(binCenters,meanValuesForBins,15^2,'r','filled')
scatter(Lprime,growthRatesMean,7^2,'k','filled')

xlim([0,40]); 
ylim([0,3]);

ylabel(['Average growth rate' 10 '(dbl/hr)']);
xlabel('Birth length (µm)');

MW_makeplotlookbetter(20);

%%
figure(hA);
%plot(Lprime,120./growthRatesMean,'r-')
plot(Lprime,(1.8)./(growthRatesMean.*log(2)),'k:','LineWidth',2);
ylim([0,60]);

saveas(hA, [OUTPUTDIR 'areaAnalysis_LvsA.svg']);
saveas(hA, [OUTPUTDIR 'areaAnalysis_LvsA.fig']);
%%

% figure; clf; hold on;
% 
% [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins]=binnedaveraging({allMotherLengths(selectedDataIdx)},{allAreasPrime(selectedDataIdx)},[0:2:40]);
% 
% scatter(allMotherLengths(selectedDataIdx),allAreasPrime(selectedDataIdx),...
%     7^2,[.7 .7 .7],'filled','MarkerFaceAlpha',.5);
% scatter(binCenters,meanValuesForBins,15^2,'k','filled')
% xlim([0,40]); %ylim([0,800]);
% 
% title('Using logarithm of length');
% ylabel('Integrated log(L(t))dt (min*µm)');
% xlabel('Birth length (µm)');
% 
% MW_makeplotlookbetter(20);

