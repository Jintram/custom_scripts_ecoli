


% A whole directory can simply be copied like this:
% copyfile 'H:\Test\FolderSource\' 'H:\Test\FolderDestination\' 

%{
% 2017
SOURCEDIR = 'H:\EXPERIMENTAL_DATA_2017\'
TARGETDIR = 'Z:\Analysis_backup\EXPERIMENTAL_DATA_2017\'

% 2016
SOURCEDIR = 'F:\EXPERIMENTAL_DATA_2014-2015_m1\'
TARGETDIR = 'Z:\Analysis_backup\EXPERIMENTAL_DATA_2014-2015_m1\'

% 2015 microscope 1
SOURCEDIR = 'F:\EXPERIMENTAL_DATA_2014-2015_m1\'
TARGETDIR = 'Z:\Analysis_backup\EXPERIMENTAL_DATA_2014-2015_m1\'

% 2015 microscope 2
SOURCEDIR = 'F:\EXPERIMENTAL_DATA_2014-2015_m2\'
TARGETDIR = 'Z:\Analysis_backup\EXPERIMENTAL_DATA_2014-2015_m2\'

%}

% Excludelist:
excludelist=[{'.','..'} arrayfun(@(x) ['pos' num2str(x)], 1:20,'UniformOutput',0) arrayfun(@(x) ['Pos' num2str(x)], 1:20,'UniformOutput',0)];

% So we need a script that takes a directory with data

% DATASET LOOP ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contentsRootDir = dir(SOURCEDIR);
sumTime=0; tic;
for idxDataset=1:numel(contentsRootDir) 
    
    % Check whether entry in root dir is folder 
    if contentsRootDir(idxDataset).isdir 
    if ~ismember(contentsRootDir(idxDataset).name,{'.','..'})    
        
        % INSIDE THE DATASET FOLDER ***************************************                
        disp(['Dataset folder found: ' contentsRootDir(idxDataset).name]);

        % Dataset dirs
        currentDatasetTargetDir = [TARGETDIR contentsRootDir(idxDataset).name '\'];
        currentDatasetSourceDir = [SOURCEDIR contentsRootDir(idxDataset).name '\'];
        % Create folder with that name in target dir            
        mkdir(currentDatasetTargetDir);
        % Get contents of dataset folder (e.g. 2017-01-06_asc990_aKG_lac)
        contentsDatasetFolder = dir([SOURCEDIR contentsRootDir(idxDataset).name]);

        % loop over sub entries of dataset (e.g. "pos1crop", "part1")
        % SEPARATE ENTRIES DATASET FOLDER==============================
        for idxSubfolder = 1:numel(contentsDatasetFolder)                        
            
            % again, check whether it is a folder that we desire
            if contentsDatasetFolder(idxSubfolder).isdir 
            if ~ismember(contentsDatasetFolder(idxSubfolder).name,excludelist)
                disp(['    Inspecting ' contentsDatasetFolder(idxSubfolder).name]);                                           
                
                % For Sub entries of subdir of dataset --------------------
                SubFolderSource = [currentDatasetSourceDir contentsDatasetFolder(idxSubfolder).name '\'];
                SubFolderTarget = [currentDatasetTargetDir contentsDatasetFolder(idxSubfolder).name '\'];
                contentsSUBSUBFolder = dir(SubFolderSource);
                mkdir(SubFolderTarget);
                
                filecount=0; foldercount=0;
                timer1start=toc;
                for SUBSUBindex = 1:numel(contentsSUBSUBFolder)
                    % Check 3rd time if folder AND if desired folder
                    if contentsSUBSUBFolder(SUBSUBindex).isdir                     
                    if ~ismember(contentsSUBSUBFolder(SUBSUBindex).name,excludelist)
                        % Copy if SUBSUB entry folder is desired
                        sourceSubSubPathToCopy = [SubFolderSource contentsSUBSUBFolder(SUBSUBindex).name];
                        targetSubSubPathToCopy = [SubFolderTarget contentsSUBSUBFolder(SUBSUBindex).name];
                        
                        % copy file
                        timer2start=toc;
                        disp(['        Copying: '     sourceSubSubPathToCopy]);
                        disp(['        Destination: ' targetSubSubPathToCopy]);
                        copyfile(sourceSubSubPathToCopy,targetSubSubPathToCopy);
                        timer2stop=toc;
                        disp(['        Done copying ' contentsSUBSUBFolder(SUBSUBindex).name '(took ' sprintf('%0.0f', (timer2stop-timer2start)/60) ' minutes).']);
                        foldercount=foldercount+1;
                        
                    else % skip if SUBSUB entry is undesired
                        if ~ismember(contentsSUBSUBFolder(SUBSUBindex).name,{'.','..'})
                            disp(['        Skipping ' SubFolderSource contentsSUBSUBFolder(SUBSUBindex).name]);        
                        end
                    end
                    else % if SUBSUB entry is file, copy
                        filecount=filecount+1;                        
                        copyfile([SubFolderSource contentsSUBSUBFolder(SUBSUBindex).name],...
                                 [SubFolderTarget contentsSUBSUBFolder(SUBSUBindex).name]);
                    end
                end                      
                % End Sub entries of subdir of dataset --------------------
            
            % some user output
            timer1stop=toc;
            disp(['    Copied ' num2str(filecount) ' files and ' num2str(foldercount) ' folders in subfolder ' contentsDatasetFolder(idxSubfolder).name]);
            disp(['    Inspection of ' contentsDatasetFolder(idxSubfolder).name ' took ' sprintf('%0.0f', (timer2stop-timer2start)/60) ' minutes']);
            
            % If it was undesired folder
            else
                if ~ismember(contentsDatasetFolder(idxSubfolder).name,{'.','..'})
                    disp(['    Skipping ' currentDatasetTargetDir contentsDatasetFolder(idxSubfolder).name]);        
                end
            end
            else
                % For dataset folder subentries, if file, just copy
                disp(['    Copying file ''' contentsDatasetFolder(idxSubfolder).name ''' to target dataset dir']);
                copyfile([currentDatasetSourceDir contentsDatasetFolder(idxSubfolder).name],...
                         [currentDatasetTargetDir contentsDatasetFolder(idxSubfolder).name]);
            end
        end % end dataset sub entry loop
        % SEPARATE ENTRIES DATASET FOLDER==================================        
        % INSIDE THE DATASET FOLDER ***************************************    
            
    end
    % For dataset entries loop:    
    else % it is a file
        % otherwise, copy file from source root to target root
        disp(['    Copying file ''' contentsRootDir(idxDataset).name ''' to target root dir']);
        copyfile([SOURCEDIR contentsRootDir(idxDataset).name],...
                 [TARGETDIR contentsRootDir(idxDataset).name]);
    end
    
end
% DATASET LOOP END+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('Yay! All done!');



