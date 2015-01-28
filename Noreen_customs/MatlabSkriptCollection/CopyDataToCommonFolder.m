% This skript reads for each position of an experiment all files from certain subfolders
% (e.g. /schnitzcells/  folder for all positions of the experiment with main folder
% /XYZ/2012-01-01/). Then all these files are copied into a new common
% folder. If wished, the filename gets a prefix with 'posX'.

% ************** ADJUST ***************************
% **NOTE: pay attention to when a backslash (\) is needed! ********

% maindir: main experimental folder: \root\date\
maindir='D:\ExperimentalDataTodo\XX\';
% useCrop: decides about which position folders to use. if 'poscrop' folders
% are to be usd: =1. if 'pos' folders: =0
useCrop=1;
% subfolder. the files of which folders should be gathered. e.g.
% \analysis\schnitzcells\
%subfolder='\images\';
% ** classical options: **
% subfolder ='\images\';
% subfolder ='\segmentation\png\';
 subfolder ='\analysis\schnitzcells\';
% subfolder ='\analysis\fluor_r\RfromFluor\';
%subfolder ='\analysis\fluor_g\g5\';

% TOCOME: specify if only a subset of files (containing a certain string) 
% from this folder should be copied into new folder

% copyOnlyLastFile: if only last file/folder of all files & folders in
% 'subfolder' should be copied: =1. Otherwise: =0. Motivation: to compare
% colony structure and size, only the lacst image in e.g. the
% /segmentation/png
% folder is relevant
copyOnlyLastFile=0;
% ** clasical options **
% colony images: =1
% analysis: =0

% addPosName: add posX to file name before copying to new folder: =1 (in
% case all files are named the same). leave name the same: =0
addPosName=1;
% ** classical choices: **
% raw images: =0
% all other stuff: =1

% combinedFolder: folder where all data has to be copied to. It's a direct
% subfolder of 'maindir'
%combinedFolder = 'AllPositions_LastmCherryImage';
% ** classical options: **
% combinedFolder = 'AllPositions_LastmCherryImage';
% combinedFolder = 'AllPositions_LastSegmentationImage';
 combinedFolder = 'AllPositions_AverageTraces';
% combinedFolder = 'AllPositions_AllFluoR6Images';
% combinedFolder = 'AllPositions_LastFluoImages';

% get all 'pos' folders (and files e.g. ('pos.STG')
allposdir=dir([maindir 'pos*']);
% sort out right directories
posfolders={};
for i=1:length(allposdir);
    currentposfolder=allposdir(i).name;
    % check if name is a directory
    checkisdir=isdir([maindir currentposfolder]);
    if checkisdir==1
        % check if 'posX' or 'posXcrop' -> search for 'crop' in name
        checkcropfolder=findstr(currentposfolder,'crop');
        if isempty(checkcropfolder) & useCrop==0
            posfolders=[posfolders, currentposfolder];
        elseif ~isempty(checkcropfolder) & useCrop==1
            posfolders=[posfolders, currentposfolder];
        end
    end
end
% -----------------------------

% make new save directory in case it does not exist yet
saveDir=[maindir, combinedFolder '\'];
if exist(saveDir)~=7
  [status,msg,id] = mkdir([saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' saveDir ' : ' msg]);
    return;
  end
end
% -----------------------------

% loop over all position folders, extract data and save into combined
% folder
for runpos=1:length(posfolders)
    disp(['Currentfolder: ' posfolders{runpos}]);
    completePathSubfolder=[maindir posfolders{runpos} subfolder];
    subdirfiles=dir(completePathSubfolder);
    % copy everything except for '.' and '..' folder
    if ~copyOnlyLastFile
        for filerun=1:length(subdirfiles)  %filerun maybe file or folder
            if strcmp(subdirfiles(filerun).name,'.')==0 & strcmp(subdirfiles(filerun).name,'..')==0
                if addPosName==1
                    newCompletePathName = [ saveDir posfolders{runpos} '_' subdirfiles(filerun).name];
                else
                    newCompletePathName = [ saveDir subdirfiles(filerun).name];
                end
            copyfile([completePathSubfolder subdirfiles(filerun).name], newCompletePathName);
            disp(['Copied file: ' subdirfiles(filerun).name])
            end
        end
    else
       filerun=length(subdirfiles);
       if filerun>0
           if strcmp(subdirfiles(filerun).name,'.')==0 & strcmp(subdirfiles(filerun).name,'..')==0
                    if addPosName==1
                        newCompletePathName = [ saveDir posfolders{runpos} '_' subdirfiles(filerun).name];
                    else
                        newCompletePathName = [ saveDir subdirfiles(filerun).name];
                    end
                copyfile([completePathSubfolder subdirfiles(filerun).name], newCompletePathName);
                disp(['Copied file: ' subdirfiles(filerun).name])
           end
       end
    end
   
end


clear i allposdir