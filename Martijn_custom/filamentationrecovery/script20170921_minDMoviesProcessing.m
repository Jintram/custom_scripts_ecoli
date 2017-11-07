

% Make sure to do both BATCH1 and BATCH2. (see commented sections)

% [6:13,14:32]

for myNum = [6:13]
%for myNum = [14:32]

    disp(['Processing ' num2str(myNum)]);
    
    TheDir = ['H:\EXPERIMENTAL_DATA_2017\2017-09-19_MinD-YFP\PART4_MinD-YFP_movies\BATCH1\m' num2str(myNum) '\'];
    %TheDir = ['H:\EXPERIMENTAL_DATA_2017\2017-09-19_MinD-YFP\PART4_MinD-YFP_movies\BATCH2\m' num2str(myNum) '\'];
    OUTDIR = 'H:\EXPERIMENTAL_DATA_2017\2017-09-19_MinD-YFP\PART4_MinD-YFP_movies\pos1\images\';

    contents = dir(TheDir);

    %%
    counter=0; mergedImage=[];
    for folderIndex=1:numel(contents)

        if any(strfind(contents(folderIndex).name,'m'))

            counter=counter+1;

            fileName = [TheDir contents(folderIndex).name];
            theImage = double(imread(fileName));

            myData(counter).theImage = theImage;

            normalizedImage = (theImage-min(theImage(:)))./(max(theImage(:))-min(theImage(:)));
            myData(counter).normalizedImage = normalizedImage;

            if counter==1
                mergedImage = normalizedImage;
            else
                mergedImage = mergedImage+normalizedImage;
            end

        end

        if any(strfind(contents(folderIndex).name,'p'))

            sourceFileName = [TheDir contents(folderIndex).name];
            copyfile(sourceFileName,[OUTDIR 'pos1-p-1-' sprintf('%03d', myNum) '.tif']); 
            disp('Copying phase');
            
        end

    end
    
    mergedImage=mergedImage./counter;

    imwrite(mergedImage, [OUTDIR 'pos1-y-' sprintf('%03d', myNum) '.tif']);
    disp('Writing yfp img');
    
end

disp('All done');
