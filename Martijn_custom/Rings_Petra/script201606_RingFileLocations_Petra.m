
% Petra's data

% settings file
load ('G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-19_part2\outputandsettings_pos1crop.mat');
% but manually overwrite some things
p.analysisDir = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-19_part2\pos1crop\analysis\';

% schnitzcells file
load ('G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-19_part2\pos1crop\data\pos1crop-Schnitz.mat');

% skeleton file
load ('G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-19_part2\pos1crop\data\pos1crop-skeletonData.mat');

% straightened fluor data file
load ('G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-19_part2\pos1crop\analysis\straightenedCells\2016-05-19pos1crop_straightFluorData.mat');
straightenedDataFolder = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-19_part2\pos1crop\analysis\straightenedCells\';
% typical filename for single cell image: 
%2016-05-19pos1crop_straightenedPlot_g_fr100cell046.tif