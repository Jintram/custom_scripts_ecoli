
% Setting _________________________________________________________________
myFolder = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\2014-10-07\468\1\';
myBFPrefix = 'BF_';
myFPPrefix = 'GFP_';

% margin to put into treshold value
tresholdMargin=0.9;

% Preparation _____________________________________________________________
myFileListing = dir(myFolder);

% Main script _____________________________________________________________

% Load file
for theFile=myFileListing'
   
    % if prefix is there
    if length(theFile.name)>=length(myBFPrefix)
    if strcmp(theFile.name(1:length(myBFPrefix)),myBFPrefix)
        % load image
        myImg=imread([myFolder theFile.name]);
        % normalize image
        myImg=double(myImg);
        myImg=(myImg-min(myImg(:)))./(max(myImg(:))-min(myImg(:)));
        % show image
        figure(1), imshow(myImg,[])
        
        % select piece for background 
        myRect = getrect();
        xmin=myRect(1); ymin=myRect(2); width=myRect(3); height=myRect(4);
        x1=xmin;y1=ymin;x2=xmin+width;y2=ymin+height;
        
        % get max background lvl        
        subImg = myImg(y1:y2,x1:x2);
        myTreshold = tresholdMargin*min(subImg(:))
        
        figure(2), imshow(subImg);
                
        % transform img
        tresholdedImg = im2bw(myImg,[],myTreshold);
        figure(3), imshow(tresholdedImg);
       
        % get edges (just to show user)
        %outlineImg = edge(tresholdedImg, 'canny');       
        [myB,myL,myN,myA] = bwboundaries(1-tresholdedImg);
        cellsAndOutlineImg = cat(3,myImg,myImg,myImg);        
        %[redCol,redRow]=find(outlineImg);        
        for compIdx=[1:myN]
        currentOutline=cell2mat(myB(compIdx));
        for idx=[1:length(currentOutline)]            
            cellsAndOutlineImg(currentOutline(idx,1),currentOutline(idx,2),:)=[1,0,0];
        end
        end
        figure(4), imshow(cellsAndOutlineImg);
        
        break; % TODO REMOVE JUST FOR TESTING        
        
    end
    end   

end

