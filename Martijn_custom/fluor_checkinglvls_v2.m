%
% Fluor level check script.
% MW 2014/10
% ===
%
% This script allows you to quantify fluor levels of bacteria in a
% microscope image. It loops over fluor files in a directory, and you need
% to select a piece which represents the background of the image. From this
% piece a treshold level is determined, and using this areas of the image 
% are identified which are fluorescent. If these areas are large enough
% they are assumed to be cells. For each of these cells the mean
% fluorescence level is determined. The signal:noise ratio is then
% determined from each cells using the area identified as background. 
%
% The Matlab script will then summarize this data for you in a plot, which
% plots all individual cells, categorized by micr. images used. Averages 
% and standard deviations are also calculated. Areas that are identified as
% cells are also highlighted in both the phase contrast and fluor image.
%
% Example input: directory with files 
% BF_30ms_2014-10-07-141019.tif
% BF_30ms_2014-10-07-141128.tif
% GFP_100ms_2014-10-07-141129.tif
% GFP_100ms_2014-10-07-141020.tif
%
% Note that images are linked together by their datenum timestamp.
%
% Output: 
% mydata.png
% mydata.xls
% mydata_mic-img-X.jpg
%
% Running the file
% ===
% Simply set parameters:
% Parameters that can be set externally
% - myFolder
% - manualLim
% - myPCPrefix and myFPPrefix
% - timeMargin
%
% And run: 
% fluor_checkinglvls_v2
%
% TODOs/note:
% - THE BACKGROUND VALUE IS DISTORTED BECAUSE I NORMALIZE THE IMAGE!
% - This script does not perform corrections to the fluor images.

disp('NOTE: Keep in mind fluor images are not corrected for lighting effects.')

% Setting _________________________________________________________________
% Make this setting such that it can be set before execution of the script.
if ~exist('myFolder','var')
    %myFolder = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\2014-10-07\negcon732\';
    %myFolder = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\2014-10-08_testicd\783\denser\';
    %myFolder = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\2014-10-07\467\1\';
    error('Please set myFolder parameter');
end
if ~exist('myPCPrefix','var') | ~exist('myFPPrefix','var')
    myPCPrefix = 'BF_';
    myFPPrefix = 'GFP_';
    disp('% myPCPrefix and myFPPrefix set to default');
end
if ~exist('timeMargin','var')
    timeMargin = 5/(24*60*60); % to determine accompanyning fluor img, 1st nr is seconds
end
ExportFilename = 'mydata';

myTextSize=8;

% margin to put into treshold value
%tresholdMargin=1.0;
if ~exist('myTresholdPercentile','var')
    myTresholdPercentile=97; % 97 for normal data, 90 for negcon
end
minCellSize=100;
limMinRatio = .9; limMaxRatio = 3.0;
limAbsoluteDifference=100;

myAnalysisFolder = [myFolder 'analysis\'];
if ~exist(myAnalysisFolder,'dir')
    mkdir(myAnalysisFolder);
end

% Preparation _____________________________________________________________
myFileListing = dir(myFolder);

multipleImgMin=[];
multipleImgMax=[];

dirMeanValues = [];

dirMeanMultipleSignalToNoise = [];
dirStdMultipleSignalToNoise = [];
dirMeanMultipleSignalMinusNoise = [];
dirStdMultipleSignalMinusNoise = [];

dirMultipleSignal={};
allMultipleSignalNoise = {};
allMultipleSignalMinusNoise = {};

if ~(myFolder(end)=='\'), disp('YOU FORGOT THE TRAILING SLASH! (in param myFolder)'), break, end

% Main script _____________________________________________________________

% Loop over files in dir 
ImagesAnalyzed=0; fluorFileNames = {};
for theFile=myFileListing'       
    
    myFileName=theFile.name;
    
    % and pick the fluor files
    %if length(myFileName)>=length(myFPPrefix) % avoid obvious non-matches
    %if strcmp(myFileName(1:length(myFPPrefix)),myFPPrefix)
    myFileName, regexptranslate('wildcard',myFPPrefix)
    regexp(myFileName, regexptranslate('wildcard',myFPPrefix))
    if ~isempty(regexp(myFileName, regexptranslate('wildcard',myFPPrefix)))
        
        % admin
        ImagesAnalyzed=ImagesAnalyzed+1;
        
        % load image
        fluorPath=[myFolder myFileName];
        myImg=imread(fluorPath);
        % normalize image
        myImg=double(myImg);
        myImgMin=min(myImg(:)); multipleImgMin(end+1)=myImgMin; % also save
        myImgMax=max(myImg(:)); multipleImgMax(end+1)=myImgMax; % also save     
        myImg=(myImg-myImgMin)./(myImgMax-myImgMin);
        % show image
        figure(1), imshow(myImg,[])
        text(10,size(myImg,2)-30,fluorPath,'Color','w','BackgroundColor','k','FontSize',myTextSize) 
        
        % select piece for background 
        disp('Select area to determine background.');
        myRect = getrect();
        
        % clear summary figure (note that getrect is the interaction w. the
        % user)
        figure(8); clf; hold on;
        
        xmin=myRect(1); ymin=myRect(2); width=myRect(3); height=myRect(4);
        x1=xmin;y1=ymin;x2=xmin+width;y2=ymin+height;        
        bgLvlImg = myImg(y1:y2,x1:x2);
        
        % select piece for treshold determination
        %{
        disp('Secect area to determine bacteria boundaries.');
        myRect = getrect();
        xmin=myRect(1); ymin=myRect(2); width=myRect(3); height=myRect(4);
        x1=xmin;y1=ymin;x2=xmin+width;y2=ymin+height;        
        tresholdLvlImg = myImg(y1:y2,x1:x2);
        %figure(2), imshow(tresholdLvlImg);
        %}
        
        % get treshold background lvl                
        %myTreshold = tresholdMargin*max(tresholdLvlImg(:))
        %myTreshold = prctile(bgLvlImg(:),myTresholdPercentile);
        % manual prctile to avoid licencing issues
        sortedbgLvlImgValues=sort(bgLvlImg(:));
        prctileidx=round(length(sortedbgLvlImgValues)*myTresholdPercentile/100);
        myTreshold=sortedbgLvlImgValues(prctileidx);
        
        % transform img
        myFilter = fspecial('average', 7)
        blurredImg = imfilter(myImg, myFilter);
        figure(99), imshow(blurredImg);        
        text(10,size(blurredImg,2)-30,fluorPath,'Color','w','BackgroundColor','k','FontSize',myTextSize) 
        tresholdedImg = im2bw(blurredImg,[],myTreshold);
        figure(3), imshow(tresholdedImg);
        text(10,size(tresholdedImg,2)-30,fluorPath,'Color','w','BackgroundColor','k','FontSize',myTextSize) 
       
        % get edges (just to show user)
        outlineImg = edge(tresholdedImg, 'canny'); % faster then bwboundaries
        %[myB,myL,myN,myA] = bwboundaries(tresholdedImg,4);
        cellsOutlineImg = cat(3,myImg,myImg,myImg);        
        [redCol,redRow]=find(outlineImg);
        for idx=[1:length(redCol)]
            cellsOutlineImg(redCol(idx),redRow(idx),:)=[1,0,0];
        end
        figure(4), imshow(cellsOutlineImg);        
        text(10,size(cellsOutlineImg,2)-30,fluorPath,'Color','w','BackgroundColor','k','FontSize',myTextSize) 
        
        % Find and show also phase contrast image to user
        % ===
        
        % find brightfield images
        % loop over files in dir
        for thePhaseFile=myFileListing'
            
            myPhaseFileName=thePhaseFile.name;

            % again, select one with prefix
            %if length(myPhaseFileName)>=length(myPCPrefix)
            %if strcmp(myPhaseFileName(1:length(myPCPrefix)),myPCPrefix)
            if ~isempty(regexp(myPhaseFileName, regexptranslate('wildcard',myPCPrefix)))
            % but now also about same time taken fluor image
            if ((theFile.datenum-thePhaseFile.datenum)<timeMargin) && ...
                    ((theFile.datenum-thePhaseFile.datenum)>-timeMargin)
                acceptedDifference=theFile.datenum-thePhaseFile.datenum
                phasePath = [myFolder myPhaseFileName];
                phaseImage = imread(phasePath);
                % normalize image
                phaseImage=double(phaseImage);
                phaseImage=(phaseImage-min(phaseImage(:)))./(max(phaseImage(:))-min(phaseImage(:)));
                % Also make on with red lines
                % ==
                % If dimensions don't match, resize outlineImg and
                % recalculate outline
                if ( length(phaseImage)/length(myImg) ) ~= 1
                    resized_outlineImg = imresize(outlineImg, length(phaseImage)/length(myImg));
                    [redCol,redRow]=find(resized_outlineImg);
                end
                % Make red border lines in phase image
                PhaseOutlineImg = cat(3,phaseImage,phaseImage,phaseImage);                        
                    for idx=[1:length(redCol)]
                        PhaseOutlineImg(redCol(idx),redRow(idx),:)=[1,0,0];
                    end
                    %figure(200), imshow(PhaseOutlineImg,[]);        
                    %text(myTextSize,size(PhaseAndOutlineImg,2)-30,fluorPath,'Color','w','BackgroundColor','k')
                break;
            end    
            end
                
        end
        figure(5), imshow(phaseImage,[]);        
        text(10,size(phaseImage,2)-30,phasePath,'Color','w','BackgroundColor','k','FontSize',myTextSize) 
        
        % Now analyze the file
        % ===
        figure(6), imshow(bgLvlImg), title('Piece of img used for mean background lvl');
        
        % determine mean of background
        meanBackground = mean(bgLvlImg(:));
        % determine mean of fluor
        fluorDataImage = double(tresholdedImg).*myImg;
        fluorValues = fluorDataImage(find(fluorDataImage>0));
        meanFluor = mean(fluorValues);
        figure(7), imshow(fluorDataImage,[]);                
               
        % get components in 
        % TODO bacterial size
        CC = bwconncomp(tresholdedImg);
        
        % make img to show data labels in
        labeledImg = cellsOutlineImg;
        figure(8); subplottight(1,2,1); 
        imshow(labeledImg,[]);        
        
        % loop over found objects (i.e. cells)
        multipleMeansFluor = []; allTxtCoordsI = []; allTxtCoordsJ = [];
        for objectIdxs = CC.PixelIdxList
            % convert to matrix
            objectIdxs=cell2mat(objectIdxs);
            % if object is large enough to be cell
            if size(objectIdxs,1)>minCellSize     
                % acquire the values at those points
                fluorValues = fluorDataImage(objectIdxs);
                % calculate mean and add to vector of means
                meanFluor = mean(fluorValues(:));
                multipleMeansFluor = [multipleMeansFluor ...
                        meanFluor];
                % get the i,j indices of the object
                [i,j] = ind2sub(size(fluorDataImage),objectIdxs);
                % text
                allTxtCoordsI=[allTxtCoordsI,i(end)];
                allTxtCoordsJ=[allTxtCoordsJ,j(end)];
                text(j(end),i(end),num2str(meanFluor),'Color','y','FontSize',myTextSize)%,'BackgroundColor','k')                 
            end
        end
        text(10,30,['meanBackground=' num2str(meanBackground)],'Color','r','FontSize',myTextSize)
        text(10,size(labeledImg,2)-30,fluorPath,'Color','w','BackgroundColor','k','FontSize',myTextSize) 
        
        % output means to user
        % ===
        % overall means
        meanBackground
        meanFluor
        signalToNoise = meanFluor/meanBackground
        
        % mean and std between cell means
        meanmultipleMeansFluor = mean(multipleMeansFluor)        
        stdmultipleMeansFluor = std(multipleMeansFluor)                
        
        % signal:noise and signal-noise for cells
        multipleSignalNoise = multipleMeansFluor./meanBackground;
        mulitipleSignalMinusNoise = multipleMeansFluor-meanBackground;
        
        % mean and std of signal:noise and signal-noise
        meanMultipleSignalToNoise = mean(multipleMeansFluor./meanBackground)
        stdMultipleSignalToNoise = std(multipleMeansFluor./meanBackground)
        meanMultipleSignalMinusNoise = mean(multipleMeansFluor-meanBackground)
        stdMultipleSignalMinusNoise = std(multipleMeansFluor-meanBackground)        
        
        % store above in arrays w. entry for each img file
        dirMultipleSignal{end+1}=multipleMeansFluor;
        dirMeanMultipleSignalToNoise(end+1) = meanMultipleSignalToNoise;
        dirStdMultipleSignalToNoise(end+1) = stdMultipleSignalToNoise;
        dirMeanMultipleSignalMinusNoise(end+1) = meanMultipleSignalMinusNoise;
        dirStdMultipleSignalMinusNoise(end+1) = stdMultipleSignalMinusNoise;
        % also store mean lvls in array
        dirMeanValues(end+1) = meanBackground;
        
        % raw data of cellular means (to make plot)
        allMultipleSignalNoise{end+1}=multipleSignalNoise;
        allMultipleSignalMinusNoise{end+1}=mulitipleSignalMinusNoise;
        
        % save filenames to recognize files later
        fluorFileNames{end+1} = myFileName;
        
        % Finish plot with summary images
        % ===
        hFig = figure(8); subplottight(1,2,2);
        imshow(PhaseOutlineImg,[]);
        text(10,size(PhaseOutlineImg,2)-30,phasePath,'Color','w','BackgroundColor','k','FontSize',myTextSize);
        %break; % TODO REMOVE JUST FOR TESTING  
        myFilePath = [myAnalysisFolder ExportFilename '_mic-img-' num2str(ImagesAnalyzed) '.jpg'];
        saveas(hFig,myFilePath)
               
    end

end

% Calculate delta value in original units
dirMeanMultipleSignalMinusNoise.*(multipleImgMax-multipleImgMin)+multipleImgMin


% Check whether something was analyzed, otherwise break analysis
if ImagesAnalyzed==0, disp('NO IMAGES FOUND TO ANALYZE!'); break; end

% Export data to Excel
myFilePath = [myAnalysisFolder ExportFilename '.xls'];
xlswrite(myFilePath,{myAnalysisFolder},'sheet1','B2');
xlswrite(myFilePath,{'Mean signal to noise'},'sheet1','B3');
xlswrite(myFilePath,dirMeanMultipleSignalToNoise,'sheet1','C3');
xlswrite(myFilePath,{'Std signal to noise'},'sheet1','B4');
xlswrite(myFilePath,dirStdMultipleSignalToNoise,'sheet1','C4');
xlswrite(myFilePath,{'Filenames'},'sheet1','B5');
xlswrite(myFilePath,fluorFileNames,'sheet1','C5');



% Make histogram of data
% ===
figure(102), [count, binLocations] = hist(myImg(:),50);
plot(binLocations,count ,'-');
title(myFolder);

% Plotting code, signal to noise
% ===

some_colors;
hFig=figure(100), clf, hold on;
spacing=.1;
for i =[1:length(allMultipleSignalNoise)]
    figure(100),plot(allMultipleSignalNoise{i},ones(length(allMultipleSignalNoise{i}),1)+spacing*(i-1),['x' mycolors(i)])
    figure(100),plot(dirMeanMultipleSignalToNoise(i),1+spacing*(i-1),['x' 'k'],'LineWidth',3)
end

plot(mean(dirMeanMultipleSignalToNoise),1-spacing,['v' 'k'],'LineWidth',3)

set(hFig, 'Units', 'pixels')

infoText=['mean over imgs=' num2str(mean(dirMeanMultipleSignalToNoise)) ' +/- ' num2str(std(dirMeanMultipleSignalToNoise)) ' (std), tresh=' num2str(myTresholdPercentile) ' (percentile)'];


ylim([1-spacing,1+spacing*(length(allMultipleSignalNoise))]);
xlim([limMinRatio,limMaxRatio])
title({myFolder,infoText,''})

myFilePath = [myAnalysisFolder ExportFilename '_signalToNoise.png'];

mySize=[4,15];
pos = get(hFig, 'Position');
set(hFig, 'Units', 'centimeters', 'Position', ...
   [2,2, ...
   mySize(2), mySize(1)]);
set(hFig, 'Units', 'centimeters', 'PaperPosition', ...
   [2,2, ...
   mySize(2), mySize(1)]);

saveas(hFig,myFilePath);

% Plotting code, signal MINUS noise
% ===

% convert back from normalized img values
if 1
    for i = [1:length(allMultipleSignalMinusNoise)]
        allMultipleSignalMinusNoise{i}=allMultipleSignalMinusNoise{i}.*(multipleImgMax(i)-multipleImgMin(i))+multipleImgMin(i);
    end
    dirMeanMultipleSignalMinusNoise=dirMeanMultipleSignalMinusNoise.*(multipleImgMax-multipleImgMin)+multipleImgMin;
end

hFig=figure(300), clf, hold on;
spacing=.1;
for i =[1:length(allMultipleSignalMinusNoise)]
    plot(allMultipleSignalMinusNoise{i},ones(length(allMultipleSignalMinusNoise{i}),1)+spacing*(i-1),['x' mycolors(i)])
    plot(dirMeanMultipleSignalMinusNoise(i),1+spacing*(i-1),['x' 'k'],'LineWidth',3)
end

plot(mean(dirMeanMultipleSignalMinusNoise),1-spacing,['v' 'k'],'LineWidth',3)
set(hFig, 'Units', 'pixels')

infoText=['mean over imgs=' num2str(mean(dirMeanMultipleSignalMinusNoise)) ' +/- ' num2str(std(dirStdMultipleSignalMinusNoise)) ' (std), tresh=' num2str(myTresholdPercentile) ' (percentile)'];

% set limits
ylim([1-spacing,1+spacing*(length(allMultipleSignalMinusNoise))]);
xlim([0,limAbsoluteDifference])

title({myFolder,infoText})

myFilePath = [myFolder 'analysis\' ExportFilename '_signalMinusNoise.png'];

mySize=[4,15];
pos = get(hFig, 'Position');
set(hFig, 'Units', 'centimeters', 'Position', ...
   [2,2, ...
   mySize(2), mySize(1)]);
set(hFig, 'Units', 'centimeters', 'PaperPosition', ...
   [2,2, ...
   mySize(2), mySize(1)]);

saveas(hFig,myFilePath);

%% EXTRA PLOT XXXXXXXXXXXXX
% ===

figure(301), clf;
set(gca,'FontSize',20);

% Administration for plotting
N = numel(dirMeanValues);
dy = 1/(N+1);
lefty = dy; righty = 1-dy;
ylocs = linspace(lefty,righty,N);

plot(dirMeanValues,ylocs,'ok','LineWidth',3);
hold on;
for i = 1:numel(dirMultipleSignal)
    plot(dirMultipleSignal{i},ones(1,numel(dirMultipleSignal{i})).*ylocs(i),'xb');
    plot(mean(dirMultipleSignal{i}),ylocs(i),'ob','LineWidth',3);
end

if exist('manualLim') % option to set ylim manually
    myLim = manualLim;
else
    myLim = max(  cell2mat( [dirMultipleSignal] )  )*1.2;
end
xlim([0,myLim]);
ylim([lefty-dy,righty+dy]);
set(gca,'YTickLabel','');

ylabel('Colonies')
xlabel('Intensity')
title(([myFolder sprintf('\n this should be done on non-normalized data!')]),'FontSize',13)

%% Build database whilst working
% ===

if ~exist('databaseValuesSignalNoise'), databaseValuesSignalNoise={}; end
databaseValuesSignalNoise{end+1}=dirMeanMultipleSignalToNoise;
if ~exist('databaseValuesMeanSignalNoise'), databaseValuesMeanSignalNoise=[]; end
databaseValuesMeanSignalNoise(end+1)=mean(dirMeanMultipleSignalToNoise);
if ~exist('databaseValuesStdSignalNoise'), databaseValuesStdSignalNoise=[]; end
databaseValuesStdSignalNoise(end+1)=std(dirMeanMultipleSignalToNoise);
if ~exist('databaseValuesNames'), databaseValuesNames={}; end
databaseValuesNames{end+1}=myFolder;


% Save whole analysis to file
save([myAnalysisFolder 'complete_analysis_' date '.mat']);



