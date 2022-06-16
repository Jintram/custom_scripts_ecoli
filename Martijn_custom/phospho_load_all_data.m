
% Phospho_load_all_data
% ===
%
% Loads all the data which is currently relevant for the regulation
% mutants.


% Set up / clear parameters
% ===
myPhosphoAuxiliary = {};
myPhosphoData = {};

% Common settings
myPhosphoAuxiliary.myRootDir = 'F:\A_Tans1_step4a_partially_analyzed_analysis\'
myPhosphoAuxiliary.myLegendNames = {}
myPhosphoAuxiliary.markers.r1 = 'o'; % this can be done prettier I guess TODO
myPhosphoAuxiliary.markers.r2 = 's';
myPhosphoAuxiliary.markers.r3 = '^';

% 732___ (Wildtype)
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos1crop', '2014-05-01','s732','r1','Wildtype','setup1','coolsnap'); 
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos4crop', '2014_06_18','s732','r2','Wildtype','setup1','coolsnap'); % 732___ colony #2 (+Tween)
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos2crop', '2014_06_18','s732','r3','Wildtype','setup1','coolsnap'); % 732___ colony #3 (+Tween)

% 733___ (pykF)
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos8crop', '2014-05-01','s733','r1','{pykF} mutant','setup1','coolsnap');
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos7crop', '2014-05-01','s733','r2','{pykF} mutant','setup1','coolsnap'); 

% 734___ (fbp)
% D:\MICROSCOPE_EXPERIMENTS\To_Analyze\2014-05-02\pos1crop\data\pos1crop-Schnitz.mat
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos1crop', '2014-05-02','s734','r1','{fbp} mutant','setup1','coolsnap');
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos2crop', '2014-05-02','s734','r2','{fbp} mutant','setup1','coolsnap');
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos4crop', '2014-05-02','s734','r3','{fbp} mutant','setup1','coolsnap');

% 735___ (pykF/ppc)
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos8crop', '2014-05-02','s735','r1','{pykF}/{ppc} mutant','setup1','coolsnap');
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos2crop', '2014_06_22','s735','r2','{pykF}/{ppc} mutant','setup1','coolsnap');
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos3crop', '2014_06_22','s735','r3','{pykF}/{ppc} mutant','setup1','coolsnap');




















