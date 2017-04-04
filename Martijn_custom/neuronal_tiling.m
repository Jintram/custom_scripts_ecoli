

%% This creates a set of 200x200 images based on a microscope dataset

theDataPath      = 'H:\EXPERIMENTAL_DATA_2017\2017-01-06_asc990_aKG_lac\';
posName          = 'pos3crop';

theImagePath = [theDataPath posName '\' 'images\'];
theSegFilePath = [theDataPath posName '\' 'segmentation\'];
dirContents = dir(theImagePath);

{dirContents.name}

filesOfInterest = find(cellfun(@(x) ~isempty(x), strfind({dirContents.name},'pos3crop-p-2-')));

frameImageName = dirContents(filesOfInterest(1)).name;
frameImagePath = [theImagePath frameImageName];
myImg = imread(frameImagePath);

% get frame number from filename
frIdxStrIndices = regexp(frameImageName,'\d');
frNrStr = frameImageName(frIdxStrIndices(3:end));

% Load segmentation file
loadedSegmentation = load([theSegFilePath posName 'seg' frNrStr '.mat'],'Lc','rect');
mySegmentation = zeros(size(myImg));
mySegmentation(loadedSegmentation.rect(1):loadedSegmentation.rect(3),loadedSegmentation.rect(2):loadedSegmentation.rect(4))=loadedSegmentation.Lc;

%% Edit segmentation, such that cell outlines get unique number and bodies get unique number
%cellBodies=zeros(size(mySegmentation));
%cellBoundaries=zeros(size(mySegmentation));
cellGradient=zeros(size(mySegmentation));
for cellNo = 1:max(mySegmentation(:))

        workImg=zeros(size(mySegmentation));
        workImg(mySegmentation==cellNo)=1;
        
        % brute force - only way I got good edges...
        h = [0 -1 0; -1 4 -1; 0 -1 0]; %h = fspecial('laplacian');  
        theBounds=imfilter(workImg,h)~=0;
        cellGradient(theBounds)=1;
        
        % Do it using imerode (doesn't work well for some reason)
        %workImg=zeros(size(mySegmentation));
        %workImg(mySegmentation==cellNo)=1;
        %cellBodies = cellBodies+imerode(workImg,strel('cube',3)); % sphere
        
        % Or do it using buondaries (also does not work well)
        %myBoundaries=bwboundaries(workImg);
        %cellBoundaries(sub2ind(size(workImg),myBoundaries{1}(:,1)',myBoundaries{1}(:,2)'))=1;
        
end

% In case using imerode (does not work well)
%mySegmentation2=mySegmentation;
%mySegmentation2(cellBodies>0)=0;
  
% Creating final image
mySegmentation2=double(mySegmentation>0);
mySegmentation2(cellGradient>0)=2;

figure; imshow(mySegmentation2,[]);

%%
figure(1); imshow(myImg,[]);
figure(2); imshow(mySegmentation,[]);

TILINGLENGTH = 200;

% Determine coordinates of tiles
xSize = size(myImg,2);
ySize = size(myImg,1);
xCoordsToSelect = [1:TILINGLENGTH:xSize];
yCoordsToSelect = [1:TILINGLENGTH:ySize];
tileCoordinateSetX = {}; tileCoordinateSetY = {};
for xidx = 1:numel(xCoordsToSelect)-1
    for yidx = 1:numel(yCoordsToSelect)-1

        tileCoordinateSetX{end+1} = [xCoordsToSelect(xidx),xCoordsToSelect(xidx+1)-1];
        tileCoordinateSetY{end+1} = [yCoordsToSelect(yidx),yCoordsToSelect(yidx+1)-1];
        
    end
end
% create extra tile for right and bottom edges (which then slightly
% overlaps)
if xSize-xCoordsToSelect(end)>TILINGLENGTH/10 & ySize-yCoordsToSelect(end)>TILINGLENGTH/10
        tileCoordinateSetX{end+1}  = [xSize-TILINGLENGTH+1,xSize];
        tileCoordinateSetY{end+1} = [ySize-TILINGLENGTH+1,ySize];
end

%%

for tileIdx = 1:numel(tileCoordinateSetX)
    tileIdx
    
    figure(3); cla; imshow(currentTile,[]);
    currentTile = myImg(tileCoordinateSetY{tileIdx}(1):tileCoordinateSetY{tileIdx}(2),tileCoordinateSetX{tileIdx}(1):tileCoordinateSetX{tileIdx}(2));

   
        
    
    pause(.1)
end

