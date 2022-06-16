

%
%THEDIR = 'H:\EXPERIMENTAL_DATA_2017\2017-09-19_MinD-YFP\PART4_MinD-YFP_movies\BATCH1\';
THEDIR = 'H:\EXPERIMENTAL_DATA_2017\2017-09-19_MinD-YFP\PART4_MinD-YFP_movies\BATCH2\';
OUTDIR = 'H:\EXPERIMENTAL_DATA_2017\2017-09-19_MinD-YFP\PART4_MinD-YFP_movies\Cropped\';


%% 
contentsMainDir = dir(THEDIR);

for mainIdx = 1:numel(contentsMainDir)
    if contentsMainDir(mainIdx).isdir & contentsMainDir(mainIdx).name(1) ~= '.' 
       
        %%
        subDirName = contentsMainDir(mainIdx).name
        
        theSubDir = [THEDIR subDirName '\'];
        contentsSubDir = dir(theSubDir);
        
        %%
        firstFlag=0;
        for subIdx = 1:numel(contentsSubDir)
            
           if contentsSubDir(subIdx).name(1) ~= '.' 
               
               %% Read the img
               theImage = imread([theSubDir contentsSubDir(subIdx).name]);
               
               % When First img in dir
               if ~firstFlag
                   firstFlag = 1;

                   %% Read img to allow crop recteangle selection                   
                   figure(1); 
                   imshow(theImage,[]);
                   theRect = int32(getrect());
                   iCrop = [theRect(2):(theRect(2)+theRect(4))];
                   jCrop = [theRect(1):(theRect(1)+theRect(3))];
                                      
                   figure(1); clf;
                   imshow(theImage(iCrop, jCrop),[]);
                   
               end                          
           
               % Crop the img
               theImageCropped = theImage(iCrop, jCrop);
               
               %% Now crop image
               myOutDir = [OUTDIR subDirName 'crop\'];
               if ~exist(myOutDir)
                   mkdir(myOutDir);
               end
               imwrite(theImageCropped,[myOutDir contentsSubDir(subIdx).name(1:end-4) '-crop.tif'])
               disp('Writing..');
               
           end

        end
        
    end
end





