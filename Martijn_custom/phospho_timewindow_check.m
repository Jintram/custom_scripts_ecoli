

% I thought there were some inconsistencies in the date, this is to check
% whether average mu determined by total length indeed fluctuates as much 
% as observed by single cell measurements.

%% Position WT (732)

p = DJK_initschnitz('pos1crop','2014-05-01','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);
% this is a bit unnessary, but quick way to get ts
[ts, ys] = ...
    DJK_plot_scatterColor(p, s_all, 'muP11_all', 'time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);

slide_plotmore(p,schnitzcells,ts)

%% Position 735

% 735___
p = DJK_initschnitz('pos8crop','2014-05-02','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);
% this is a bit unnessary, but quick way to get ts
[ts, ys] = ...
    DJK_plot_scatterColor(p, s_all, 'muP11_all', 'time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 1);

slide_plotmore(p,schnitzcells,ts)