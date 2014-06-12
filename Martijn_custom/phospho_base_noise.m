

% Script that collects all necessary data from saved Schnitzes, and then
% uses a few functions to process this data and make plots.
%
% These plots aim to quantify the noise of the cells.
%
% MW 2014/06/12
%
% Previous version: phospho_scripts_temp_mutants_obtaining_scatter.m


%% Settings

some_colors;

% Common settings
myPhosphoData.myRootDir = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\'
myPhosphoData.myLegendNames = {}

%{
for groupname={'s732','s733','s734','s735','s736'}
    groupname=num2str(cell2mat(groupname))
    zz.(groupname) = {}
end
%}

%% Load data ______________________________________________________________



% 732___
myPhosphoData=phospho_loadscatterdata(myPhosphoData, 'pos1crop', '2014-05-01','s732','r1','Wildtype');

% 733___
myPhosphoData=phospho_loadscatterdata(myPhosphoData, 'pos8crop', '2014-05-01','s733','r1','\Delta{pykF} mutant');

% 735___
myPhosphoData=phospho_loadscatterdata(myPhosphoData, 'pos8crop', '2014-05-02','s735','r1','\Delta{pykF}/\Delta{ppc} mutant');


%% Plot data ______________________________________________________________

nrpositions = 3;
mymap = colorGray(nrpositions);
fignumstart=1;
figure(fignumstart); clf; figure(fignumstart+1); clf; figure(fignumstart+2); clf;
lines1=[];lines2=[];lines3=[]; % stores line handles for legends

% Plot all 
[lgh1,lgh2,lgh3] = sliding_window_stds_plot(myPhosphoData.('s732').('r1').xvalues,myPhosphoData.('s732').('r1').yvalues,fignumstart,mymap(1,:),1,1);
lines1 = [lines1 lgh1]; lines2 = [lines2 lgh2]; lines3 = [lines3 lgh3];
[lgh1,lgh2,lgh3] = sliding_window_stds_plot(myPhosphoData.('s733').('r1').xvalues,myPhosphoData.('s733').('r1').yvalues,fignumstart,mymap(2,:),1,1);
lines1 = [lines1 lgh1]; lines2 = [lines2 lgh2]; lines3 = [lines3 lgh3];
[lgh1,lgh2,lgh3] = sliding_window_stds_plot(myPhosphoData.('s735').('r1').xvalues,myPhosphoData.('s735').('r1').yvalues,fignumstart,mymap(3,:),1,1);
lines1 = [lines1 lgh1]; lines2 = [lines2 lgh2]; lines3 = [lines3 lgh3];

% add legend
figure(fignumstart); legend(lines1, myPhosphoData.myLegendNames,'Location','best');
figure(fignumstart+1); legend(lines2, myPhosphoData.myLegendNames,'Location','best');
figure(fignumstart+2); legend(lines3, myPhosphoData.myLegendNames,'Location','best');

















