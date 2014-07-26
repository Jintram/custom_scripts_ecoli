

% Loads all the data which is currently relevant for the regulation
% mutants.



% 732___
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos1crop', '2014-05-01','s732','r1','Wildtype'); 
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos4crop', '2014_06_18','s732','r2','Wildtype'); % 732___ colony #2 (+Tween)
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos2crop', '2014_06_18','s732','r3','Wildtype'); % 732___ colony #3 (+Tween)

% 733___
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos8crop', '2014-05-01','s733','r1','{pykF} mutant');

% 734___
% D:\MICROSCOPE_EXPERIMENTS\To_Analyze\2014-05-02\pos1crop\data\pos1crop-Schnitz.mat
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos1crop', '2014-05-02','s734','r1','{fbp} mutant');
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos2crop', '2014-05-02','s734','r2','{fbp} mutant');
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos4crop', '2014-05-02','s734','r3','{fbp} mutant');

% 735___
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos8crop', '2014-05-02','s735','r1','{pykF}/{ppc} mutant');
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos2crop', '2014_06_22','s735','r2','{pykF}/{ppc} mutant');
[myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData, myPhosphoAuxiliary, 'pos3crop', '2014_06_22','s735','r3','{pykF}/{ppc} mutant');



