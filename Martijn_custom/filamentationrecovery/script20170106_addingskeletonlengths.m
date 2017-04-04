


% This script adds the skeleton lengths to the Rutger data

% Currently the algorithm spits out many warning because extrapolation is
% going wrong. Postponed this procedure for now.

% General parameter
p.micronsPerPixel =     0.04065; % old camera (coolsnap)

%% Dataset 2013-12-09\pos3crop

% specific to data
BASEDIR = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-12-09_Rutger2_TET\';
p.movieName =           'pos3crop';
p.tracksDir =           [BASEDIR p.movieName '\data\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\data\';
p.analysisDir =         [BASEDIR p.movieName '\analysis\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\analysis\';
p.segmentationDir =     [BASEDIR p.movieName '\segmentation\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\segmentation\';

% Fix frame numbers
MW_addframenrsfix
% Add the lengths
NDL_addToSchnitzes_skeletonLengthMW(p); % saves schnitzcells param

%% For dataset 2013-09-24\pos4crop

% specific to data
BASEDIR = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\';
p.movieName =           'pos4crop';
p.tracksDir =           [BASEDIR p.movieName '\data\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\data\';
p.analysisDir =         [BASEDIR p.movieName '\analysis\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\analysis\';
p.segmentationDir =     [BASEDIR p.movieName '\segmentation\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\segmentation\';

% Fix frame numbers
MW_addframenrsfix
% Add the lengths
NDL_addToSchnitzes_skeletonLengthMW(p); % saves schnitzcells param

%% Dataset \2013-12-16\pos4crop

% specific to data
BASEDIR = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-12-16_Rutger3_TET\';
p.movieName =           'pos4crop';
p.tracksDir =           [BASEDIR p.movieName '\data\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\data\';
p.analysisDir =         [BASEDIR p.movieName '\analysis\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\analysis\';
p.segmentationDir =     [BASEDIR p.movieName '\segmentation\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\segmentation\';

% Fix frame numbers
MW_addframenrsfix
% Add the lengths
NDL_addToSchnitzes_skeletonLengthMW(p); % saves schnitzcells param

%% Dataset \2013-09-24\pos5crop

% specific to data
BASEDIR = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\';
p.movieName =           'pos5crop';
p.tracksDir =           [BASEDIR p.movieName '\data\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\data\';
p.analysisDir =         [BASEDIR p.movieName '\analysis\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\analysis\';
p.segmentationDir =     [BASEDIR p.movieName '\segmentation\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\segmentation\';

% Fix frame numbers
MW_addframenrsfix
% Add the lengths
NDL_addToSchnitzes_skeletonLengthMW(p); % saves schnitzcells param

%% Dataset \2013-12-16\pos5crop

% specific to data
BASEDIR = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-12-16_Rutger3_TET\';
p.movieName =           'pos5crop';
p.tracksDir =           [BASEDIR p.movieName '\data\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\data\';
p.analysisDir =         [BASEDIR p.movieName '\analysis\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\analysis\';
p.segmentationDir =     [BASEDIR p.movieName '\segmentation\']; 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\segmentation\';

% Fix frame numbers
MW_addframenrsfix
% Add the lengths
NDL_addToSchnitzes_skeletonLengthMW(p); % saves schnitzcells param

%%
%{
        ... Note that the parameter ONESTOANALYZE makes a subselection of this data.
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos3_long.mat',...             1
        ... ^ Comes from ..\2013-12-09\pos3crop\data\
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos4.mat',...                  2
        ... ^ Comes from ..\2013-09-24\pos4crop\data\ (ugly images)
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos4_long.mat',...             3
        ... ^ Comes from ..\2013-12-16\pos4crop\data\ (ugly images)
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos5.mat',...                  4
        ... ^ Comes from ..\2013-09-24\pos5crop\data\ (very little data)
        'G:\FilamentationRecoveryData\Rutger\F schijf AmolfBackup_3april2014\USE_DIV\1uM_pos5_long.mat',...             5
        ... ^ Comes from ..\2013-12-16\pos5crop\data\ (ugly images)
%}