
% Setting _________________________________________________________________
myFolder = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\2014-10-07\453\1\';
myPCPrefix = 'BF_';
myFPPrefix = 'GFP_';
ExportFilename = 'mydata';

% margin to put into treshold value
%tresholdMargin=1.0;
myTresholdPercentile=97;
minCellSize=100;
timeMargin = 5/(24*60*60); % to determine accompanyning fluor img
limMinRatio = .9; limMaxRatio = 3.0;

% Preparation _____________________________________________________________
myFileListing = dir(myFolder);

dirMeanMultipleSignalToNoise = [];
dirStdMultipleSignalToNoise = [];
allMultipleSignalNoise = {};

% Main script _____________________________________________________________

% Loop over files in dir 
for theFile=myFileListing'
   
    myFileName=theFile.name;
    
    % and pick the fluor files
    if length(myFileName)>=length(myFPPrefix) % avoid obvious non-matches
    if strcmp(myFileName(1:length(myFPPrefix)),myFPPrefix)
        % load image
        fluorPath=[myFolder myFileName];
        myImg=imread(fluorPath);
        % normalize image
        myImg=double(myImg);
        myImg=(myImg-min(myImg(:)))./(max(myImg(:))-min(myImg(:)));
        % show image
        figure(1), imshow(myImg,[])
        text(10,size(phaseImage,2)-30,fluorPath,'Color','w','BackgroundColor','k') 
        
        % select piece for background 
        disp('Select area to determine background.');
        myRect = getrect();
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
        myTreshold = prctile(tresholdLvlImg(:),myTresholdPercentile);
                
        % transform img
        myFilter = fspecial('average', 7)
        blurredImg = imfilter(myImg, myFilter);
        figure(99), imshow(blurredImg);        
        text(10,size(phaseImage,2)-30,fluorPath,'Color','w','BackgroundColor','k') 
        tresholdedImg = im2bw(blurredImg,[],myTreshold);
        figure(3), imshow(tresholdedImg);
        text(10,size(phaseImage,2)-30,fluorPath,'Color','w','BackgroundColor','k') 
       
        % get edges (just to show user)
        outlineImg = edge(tresholdedImg, 'canny'); % faster then bwboundaries
        %[myB,myL,myN,myA] = bwboundaries(tresholdedImg,4);
        cellsAndOutlineImg = cat(3,myImg,myImg,myImg);        
        [redCol,redRow]=find(outlineImg);
        for idx=[1:length(redCol)]
            cellsAndOutlineImg(redCol(idx),redRow(idx),:)=[1,0,0];
        end
        figure(4), imshow(cellsAndOutlineImg);        
        text(10,size(phaseImage,2)-30,fluorPath,'Color','w','BackgroundColor','k') 
        
        % Find and show also phase contrast image to user
        % ===
        
        % find brightfield image
        % loop over files in dir
        for thePhaseFile=myFileListing'
            
            myPhaseFileName=thePhaseFile.name;

            % again, select one with prefix
            if length(myPhaseFileName)>=length(myPCPrefix)
            if strcmp(myPhaseFileName(1:length(myPCPrefix)),myPCPrefix)
            % but now also about same time taken fluor image
            if(thePhaseFile.datenum-theFile.datenum)<timeMargin
                phasePath = [myFolder myPhaseFileName];
                phaseImage = imread(phasePath);
                break;
            end    
            end
            end
                
        end
        figure(5), imshow(phaseImage,[]);        
        text(10,size(phaseImage,2)-30,phasePath,'Color','w','BackgroundColor','k') 
        
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
        labeledImg = cellsAndOutlineImg;
        figure(8), imshow(labeledImg,[]);
        
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
                text(j(end),i(end),num2str(meanFluor),'Color','y')%,'BackgroundColor','k')                 
            end
        end

        % output means to user
        meanBackground
        meanFluor
        signalToNoise = meanFluor/meanBackground
        
        meanmultipleMeansFluor = mean(multipleMeansFluor)        
        stdmultipleMeansFluor = std(multipleMeansFluor)
        
        MultipleSignalNoise = multipleMeansFluor./meanBackground;
        
        meanMultipleSignalToNoise = mean(multipleMeansFluor./meanBackground)
        stdMultipleSignalToNoise = std(multipleMeansFluor./meanBackground)
        
        dirMeanMultipleSignalToNoise=[dirMeanMultipleSignalToNoise,meanMultipleSignalToNoise];
        dirStdMultipleSignalToNoise=[dirStdMultipleSignalToNoise,stdMultipleSignalToNoise];
        
        allMultipleSignalNoise{end+1}=MultipleSignalNoise;
        
        %break; % TODO REMOVE JUST FOR TESTING        
               
    end
    end   

end

% Export data to Excel
myFilePath = [myFolder ExportFilename '.xls']
xlswrite(myFilePath,{myFolder},'sheet1','B2');
xlswrite(myFilePath,{'Mean signal to noise'},'sheet1','B3');
xlswrite(myFilePath,dirMeanMultipleSignalToNoise,'sheet1','C3');
xlswrite(myFilePath,{'Std signal to noise'},'sheet1','B4');
xlswrite(myFilePath,dirStdMultipleSignalToNoise,'sheet1','C4');



some_colors;
hFig=figure(100), clf, hold on;
spacing=.1;
for i =[1:length(allMultipleSignalNoise)]
    figure(100),plot(allMultipleSignalNoise{i},ones(length(allMultipleSignalNoise{i}),1)+spacing*(i-1),['x' mycolors(i)])
end
ylim([1-spacing,1+spacing*(length(allMultipleSignalNoise))]);
xlim([limMinRatio,limMaxRatio])
title(myFolder)

myFigFilePath = [myFolder ExportFilename '.png'];

size=[100,500];
pos = get(hFig, 'Position');
set(hFig, 'Units', 'pixels', 'Position', ...
   [pos(1), pos(2)+pos(4)-size(1), ...
   size(2), size(1)]);
set(hFig, 'Units', 'pixels', 'PaperPosition', ...
   [pos(1), pos(2)+pos(4)-size(1), ...
   size(2), size(1)]);

saveas(hFig,myFigFilePath);







