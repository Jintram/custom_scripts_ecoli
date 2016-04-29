



%%
PEAKTRESHOLD=400;
figure(1); clf;
figure(101); clf;
plotcolors = 'rrr';
scatterX={}; scatterY={}; 


%%
for datasetIndex = 2:3

    %%
    %load ([p.tracksDir p.movieName '-skeletonData.mat']);
    load (['G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos' num2str(datasetIndex) 'crop\data\pos' num2str(datasetIndex) 'crop-skeletonData.mat']);

    %%
    %load(saveLocationMatFile([p.analysisDir 'straightenedPlots\' p.movieDate p.movieName '_straightFluorData.mat']);
    load (['G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos' num2str(datasetIndex) 'crop\analysis\straightenedCells\2016-04-07pos' num2str(datasetIndex) 'crop_straightFluorData.mat']);


    %%

    figure(1); hold on;
    %clf;
    %hold on;

    pileofXmicron=[];
    pileofXPx = [];
    pileofmeanY=[];
    selectionVector=[];
    pileofLengthsOfBacteriaInMicrons=[];
    correspondingLengthsToLocations=[];
    correspondingLengthsToLocationsPx=[];
    for framenr=1:numel(allpeakXMicrons)

        meanYFrame  = allpeakmeanY{framenr};
        xpeaksFrame  = allpeakXMicrons{framenr};
        lengthsFrame = allLengthsOfBacteriaInMicrons{framenr};

        xpeaksFramePx  = allpeakXPixels{framenr};
        lengthsFramePx = allLengthsOfBacteriaInPixels{framenr};
        
        %disp(['adding ' num2str(numel(xpeaksFrame)) ' cells for frame ' num2str(framenr) '..']);
        pileofLengthsOfBacteriaInMicrons =...
            [pileofLengthsOfBacteriaInMicrons lengthsFrame];

        for cellno = 1:numel(xpeaksFrame)

            %plot(xpeaksFrame)
            if ~isempty(xpeaksFrame)

                % raw data
                pileofXmicron = [pileofXmicron xpeaksFrame{cellno}'];
                pileofXPx     = [pileofXPx xpeaksFramePx{cellno}'];
                fluorPeaksThisCell = meanYFrame{cellno};
                pileofmeanY = [pileofmeanY fluorPeaksThisCell];

                % create duplicate lengths for datapoints belonging to same
                % bacteria
                duplicatedLengths = ones(1,numel(xpeaksFrame{cellno}))*lengthsFrame(cellno);
                correspondingLengthsToLocations = [correspondingLengthsToLocations duplicatedLengths];
                % same w. pixels
                duplicatedLengthsPx = ones(1,numel(xpeaksFrame{cellno}))*lengthsFramePx(cellno);
                correspondingLengthsToLocationsPx = [correspondingLengthsToLocationsPx duplicatedLengthsPx];
                
                % select peaks that are >200
                selectionVector = [selectionVector fluorPeaksThisCell>PEAKTRESHOLD];

            end
            %selectedataXmicron = ...

        end

    end

    %plot(correspondingLengthsToLocations,pileofXmicron,'.')
    %plot(correspondingLengthsToLocations,pileofXmicron./correspondingLengthsToLocations,'.')
    indicesToSelect = find(selectionVector);
    scatterX{datasetIndex}=correspondingLengthsToLocations(indicesToSelect);
    scatterY{datasetIndex}=pileofXmicron(indicesToSelect)./correspondingLengthsToLocations(indicesToSelect);
    plot(scatterX{datasetIndex},scatterY{datasetIndex},'x','Color',plotcolors(datasetIndex))

    % and in Pixels
    figure(101); hold on;
    scatterXPx{datasetIndex}=correspondingLengthsToLocationsPx(indicesToSelect);
    scatterYPx{datasetIndex}=pileofXPx(indicesToSelect)./correspondingLengthsToLocationsPx(indicesToSelect);
    plot(scatterXPx{datasetIndex},scatterYPx{datasetIndex},'.','Color',plotcolors(datasetIndex))
    
    
    %plot(allpeak)

end

%% Sanity check
MICRONSPERPIXEL=0.0431;

figure(1); hold on;
plot(scatterXPx{datasetIndex}*MICRONSPERPIXEL,scatterYPx{datasetIndex},'.b')

%% 
NRCONTOURLINES=3;
MAKEUSEOFSYMMETRY=1;
COLOR='b';

figure(2); clf; hold on;

if ~MAKEUSEOFSYMMETRY
    data = [[scatterX{:}]; [scatterY{:}]]';
else
    data = [[[scatterX{:}] [scatterX{:}]]; [[scatterY{:}], 1-[scatterY{:}]]]'; % symmetry data
end

plot([scatterX{:}],[scatterY{:}],['.' COLOR],'MarkerSize',7);
if MAKEUSEOFSYMMETRY
    plot([scatterX{:}],1-[scatterY{:}],['.' COLOR],'MarkerSize',7); % plot symmetric
end
[bandwidth,density,X,Y] = kde2d(data);      
%[C, l1] = contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);

% density corrected for data loss at higher lengths..
%[C, l1] = contour(X,Y,density.*X,NRCONTOURLINES,'-k','LineWidth',2);
[C, l1] = contour(X,Y,density.*X,NRCONTOURLINES,'-k','LineWidth',2);

%title('kde, density directed (pdf(x)*x)')
xlabel(['Length of cell [' '\mu' 'm]']);
ylabel('ftsA peaks location');
%figure(); hist([scatterX{:}])

maxX = max([scatterX{:}]);
xlim([0, maxX]);
ylim([0,1]);

MW_makeplotlookbetter(20);

% plot helping lines at 1/2n
N=5;
for i=1:N
    for j = 1:(i*2-1)
        plot([0, maxX], [(j)/(2*i) (j)/(2*i)],'-','Color',[.5 .5 .5],'LineWidth',N-i+1)
    end
end

%% Create overlay plot if division ratios are avaible
NEWMAXX = 50;

figure(3); clf; hold on;
[C, l1] = contour(X,Y,density.*X,NRCONTOURLINES,'-k','LineWidth',2);

% data that should be obtained from
% script20160429_filamentRecoveryDivisionRatioss
if exist('Ratios','var')
    for datasetIdx = 1:numel(datasetsPaths)
        plot(myLengthSumNewborns{datasetIdx},Ratios{datasetIdx},'o', 'Color', PLOTCOLORS(datasetIdx,:),'LineWidth',2);
    end
end

xlim([0, NEWMAXX]);
ylim([0,1]);

xlabel(['Length of cell [' '\mu' 'm]']);
ylabel('ftsA peaks location / division location');

MW_makeplotlookbetter(15)

%% 

figure(4);clf;

[n,c] = hist(pileofLengthsOfBacteriaInMicrons,50);
plot(c,n./(sum(n)*(c(2)-c(1))),'o-k','MarkerSize',15,'LineWidth',3);
MW_makeplotlookbetter(20);
xlabel('Cell length [um]');
ylabel('Probability');

%% more sanity checks

%load (['G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos' num2str(datasetIndex) 'crop\data\pos' num2str(datasetIndex) 'crop-Schnitz.mat']);







