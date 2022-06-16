
% Setting _________________________________________________________________
myFolder = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\2014-10-07\453\1\';
myPCPrefix = 'BF_';
myFPPrefix = 'GFP_';

% margin to put into treshold value
tresholdMargin=1.0;
myTresholdPercentile=97;
timeMargin = 5/(24*60*60); % to determine accompanyning fluor img

% Preparation _____________________________________________________________
myFileListing = dir(myFolder);

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
        disp('Secect area to determine bacteria boundaries.');
        myRect = getrect();
        xmin=myRect(1); ymin=myRect(2); width=myRect(3); height=myRect(4);
        x1=xmin;y1=ymin;x2=xmin+width;y2=ymin+height;        
        tresholdLvlImg = myImg(y1:y2,x1:x2);        
        %figure(2), imshow(tresholdLvlImg);
        
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
        
        break; % TODO REMOVE JUST FOR TESTING        
        
    end
    end   

end

