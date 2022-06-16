

%% Initialize experimental details, dataset 1
%{
baseDir = 'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-09\pos3crop\';
p.movieName = 'pos3crop';
p.movieDate = '2013-12-09';
p.micronsPerPixel= 0.04065;
p.tracksDir = [baseDir 'data\'];
p.segmentationDir= [baseDir 'segmentation\'];
p.analysisDir= [baseDir 'analysis\'];
%}

%% Initialize experimental details, dataset 2
% REDONE
%{
baseDir = 'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos4crop\';
p.movieName = 'pos4crop';
p.movieDate = '2013-09-24';
p.micronsPerPixel= 0.04065;
p.tracksDir = [baseDir 'data\'];
p.segmentationDir= [baseDir 'segmentation\'];
p.analysisDir= [baseDir 'analysis\'];
%}

%% Initialize experimental details, dataset 3
% REDONE
%{
baseDir = 'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos4crop\';
p.movieName = 'pos4crop';
p.movieDate = '2013-12-16';
p.micronsPerPixel= 0.04065;
p.tracksDir = [baseDir 'data\'];
p.segmentationDir= [baseDir 'segmentation\'];
p.analysisDir= [baseDir 'analysis\'];
%}

%% Initialize experimental details, dataset 4
%{
baseDir = 'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-09-24\pos5crop\';
p.movieName = 'pos5crop';
p.movieDate = '2013-09-24';
p.micronsPerPixel= 0.04065;
p.tracksDir = [baseDir 'data\'];
p.segmentationDir= [baseDir 'segmentation\'];
p.analysisDir= [baseDir 'analysis\'];
%}

%% Initialize experimental details, dataset 5
% I think this was also redone
%{
% When performing this analysis, I noticed the following..
% This dataset has some segmentation issues: there is one cell that seemed
% to have died that is segmented inproperly, and also there are artefacts
% in the last frames that are segmentented. Importantly, there in 
% e.g. frame 1179 there is a "cut line" running through the whole colony,
% introducing artefacts in large numbers of cells. Start of dataset with
% actual divisions of filaments seems to be OK on first glance though.\
%
% Therefor, for this dataset, I manually executed 
% NDL_addToSchnitzes_skeletonLengthMW with an adjusted framerange 
% only up to 1177, and also manually executed lines below to achieve this
% frameRange.
% 
% There were also some frames missing, so I used
% myfrms=unique(schnitzcells);
% frameRange=myfrms(1:find(myfrms==1177));
% -MW

baseDir = 'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\2013-12-16\pos5crop\';
p.movieName = 'pos5crop';
p.movieDate = '2013-12-16';
p.micronsPerPixel= 0.04065;
p.tracksDir = [baseDir 'data\'];
p.segmentationDir= [baseDir 'segmentation\'];
p.analysisDir= [baseDir 'analysis\'];
%}

%% Update schnitzcells with frame_nrs if necessary
load([p.tracksDir p.movieName '-Schnitz.mat']);
schnitzcells = MW_calculateframe_nrs(schnitzcells);
save([p.tracksDir p.movieName '-Schnitz.mat'],'schnitzcells');

%% Start analysis
NDL_addToSchnitzes_skeletonLengthMW(p);

frameRange = unique([schnitzcells.frame_nrs]);
fluorColor='g';
MW_straightenbacteria(p, frameRange, fluorColor);




