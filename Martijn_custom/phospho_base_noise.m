

% Script that collects all necessary data from saved Schnitzes, and then
% uses a few functions to process this data and make plots.
%
% These plots aim to quantify the noise of the cells.
%
% MW 2014/06/12
%
% Previous version: phospho_scripts_temp_mutants_obtaining_scatter.m


%% Settings and loading of data.

% Load color scheme
some_colors;

%{
for groupname={'s732','s733','s734','s735','s736'}
    groupname=num2str(cell2mat(groupname))
    zz.(groupname) = {}
end
%}

% Loading of data
% ===
phospho_load_all_data % for convenience: edit phospho_load_all_data


%% Plot mean, std.dev. and noise over time ________________________________
% Per colony
divbymean = 1;

nrpositions = 3;
mymap = colorGray(nrpositions);
fignumstart=1;
figure(fignumstart); clf; figure(fignumstart+1); clf; figure(fignumstart+2); clf; figure(fignumstart+3); clf;
lines1=[];lines2=[];lines3=[]; lines4=[];% stores line handles for legends

% Plot all 
% 732___
% (Different plots are made in the sliding_window_stds_plot function.)
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s732').('r1').xvalues,myPhosphoData.('s732').('r1').yvalues,fignumstart,somemarkers(1),preferredcolors(1,:),1,divbymean,'r1');
lines1 = [lines1 lgh1]; lines2 = [lines2 lgh2]; lines3 = [lines3 lgh3]; lines4 = [lines4 lgh4]; % add first line to list
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s732').('r2').xvalues,myPhosphoData.('s732').('r2').yvalues,fignumstart,somemarkers(1),preferredcolors(1,:),1,divbymean,'r2');
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s732').('r3').xvalues,myPhosphoData.('s732').('r3').yvalues,fignumstart,somemarkers(1),preferredcolors(1,:),1,divbymean,'r3');

% Data on other mutant (733)
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s733').('r1').xvalues,myPhosphoData.('s733').('r1').yvalues,fignumstart,somemarkers(2),preferredcolors(2,:),1,divbymean,'r1');
lines1 = [lines1 lgh1]; lines2 = [lines2 lgh2]; lines3 = [lines3 lgh3]; lines4 = [lines4 lgh4]; % add first line to list
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s733').('r2').xvalues,myPhosphoData.('s733').('r2').yvalues,fignumstart,somemarkers(2),preferredcolors(2,:),1,divbymean,'r2');

% 734___
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s734').('r1').xvalues,myPhosphoData.('s734').('r1').yvalues,fignumstart,somemarkers(3),preferredcolors(3,:),1,divbymean,'r1');
lines1 = [lines1 lgh1]; lines2 = [lines2 lgh2]; lines3 = [lines3 lgh3]; lines4 = [lines4 lgh4]; % add first line to list
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s734').('r2').xvalues,myPhosphoData.('s734').('r2').yvalues,fignumstart,somemarkers(3),preferredcolors(3,:),1,divbymean,'r2');
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s734').('r3').xvalues,myPhosphoData.('s734').('r3').yvalues,fignumstart,somemarkers(3),preferredcolors(3,:),1,divbymean,'r3');

% 735___
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s735').('r1').xvalues,myPhosphoData.('s735').('r1').yvalues,fignumstart,somemarkers(4),preferredcolors(4,:),1,divbymean,'r1');
lines1 = [lines1 lgh1]; lines2 = [lines2 lgh2]; lines3 = [lines3 lgh3]; lines4 = [lines4 lgh4]; % add first line to list
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s735').('r2').xvalues,myPhosphoData.('s735').('r2').yvalues,fignumstart,somemarkers(4),preferredcolors(4,:),1,divbymean,'r2');
[lgh1,lgh2,lgh3,lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,myPhosphoData.('s735').('r3').xvalues,myPhosphoData.('s735').('r3').yvalues,fignumstart,somemarkers(4),preferredcolors(4,:),1,divbymean,'r3');


% add legend
theLegendNames = unique(myPhosphoAuxiliary.myLegendNames,'stable'); % b/c repetitions not in list of lines, get unique names.
figure(fignumstart); legend(lines1, theLegendNames,'Location','best');
figure(fignumstart+1); legend(lines2, theLegendNames,'Location','best');
figure(fignumstart+2); legend(lines3, theLegendNames,'Location','best');
figure(fignumstart+3); legend(lines4, theLegendNames,'Location','NorthEast');

% Some markup for "article" writing course
figure(fignumstart+3), xlim([0.3 1.2]), ylim([0 .7]), 

%% Create (mu,noise) plot with Schnitzes as datapoints ____________________

% Process data

% TODO: should be updated such that looping goes automatically, and
% phospho_mu_vs_noise shouldn't need to require anything else but
% myPhosphoData, myPhosphoAuxiliary and strain and repetition.

% 732 data (,'Wildtype')
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s732','r1');
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s732','r2');
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s732','r3');

% 733 data (pykF)
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s733','r1');
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s733','r2');

% 734 data (fbp)
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s734','r1');
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s734','r2');
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s734','r3');

% 735 data (pykF/ppc)
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s735','r1');
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s735','r2');
myPhosphoData = phospho_mu_vs_noise(myPhosphoData,'s735','r3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that here we need to inspect trend INSIDE colonies, so we cannot
% group colonies together!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create fit lines
x_values = [0:0.2:2];
% 732
std_fitCoef1_732 = polyfit(myPhosphoData.('s732').('r1').the_means,myPhosphoData.('s732').('r1').the_noises,1);
std_fitted_732 = std_fitCoef1_732(1)*x_values + std_fitCoef1_732(2);
% 733
std_fitCoef1_733 = polyfit(myPhosphoData.('s733').('r1').the_means,myPhosphoData.('s733').('r1').the_noises,1);
std_fitted_733 = std_fitCoef1_733(1)*x_values + std_fitCoef1_733(2);
% 734
std_fitCoef1_734 = polyfit(myPhosphoData.('s734').('r1').the_means,myPhosphoData.('s734').('r1').the_noises,1);
std_fitted_734 = std_fitCoef1_734(1)*x_values + std_fitCoef1_734(2);
% 735
std_fitCoef1_735 = polyfit(myPhosphoData.('s735').('r1').the_means,myPhosphoData.('s735').('r1').the_noises,1);
std_fitted_735 = std_fitCoef1_735(1)*x_values + std_fitCoef1_735(2);
% Overall fit
%{
all_means = [all_means_732,all_means_734,all_means_735];
all_noises = [all_noises_732,all_noises_734,all_noises_735];
std_fitCoef1_all = polyfit(all_means,all_noises,1);
std_fitted_all = std_fitCoef1_all(1)*x_values + std_fitCoef1_all(2);
%}


% Plot data

nrpositions = 3;
mymap = colorGray(nrpositions);
h = figure(5); clf; hold on;
set(h, 'Position', [100, 100, 800+100, 600+100])
set(gca,'FontSize',20);
%title('Growth speed vs. noise (1 point = 1 ''individual'')');
xlabel('\mu (doublings/hr)');
ylabel('Noise'); % == std. dev. / mean
%axis([0 2 0 2])
lines4=[];

% Plot all 
% 732
plot(myPhosphoData.('s732').('r1').the_means,myPhosphoData.('s732').('r1').the_noises,'.','color',preferredcolors(1,:),'MarkerSize',7)
lgh4 = plot(x_values,std_fitted_732,'--','LineWidth',4,'color',preferredcolors(1,:),'Marker',somemarkers(1)); % fit line
lines4 = [lines4 lgh4];  % legend line
% 733
plot(myPhosphoData.('s733').('r1').the_means,myPhosphoData.('s733').('r1').the_noises,'.','color',preferredcolors(2,:),'MarkerSize',7)
lgh4 = plot(x_values,std_fitted_733,'--','LineWidth',4,'color',preferredcolors(2,:),'Marker',somemarkers(2)); % fit line
lines4 = [lines4 lgh4]; % legend line
% 734
plot(myPhosphoData.('s734').('r1').the_means,myPhosphoData.('s734').('r1').the_noises,'.','color',preferredcolors(3,:),'MarkerSize',7)
lgh4 = plot(x_values,std_fitted_734,'--','LineWidth',4,'color',preferredcolors(3,:),'Marker',somemarkers(3)); % fit line
lines4 = [lines4 lgh4]; % legend line
% 735
plot(myPhosphoData.('s735').('r1').the_means,myPhosphoData.('s735').('r1').the_noises,'.','color',preferredcolors(4,:),'MarkerSize',7)
lgh4 = plot(x_values,std_fitted_735,'--','LineWidth',4,'color',preferredcolors(4,:),'Marker',somemarkers(4)); % fit line
lines4 = [lines4 lgh4]; % legend line

% edit limits 
%[all_means_732 all_means_733 all_means_734 all_means_735]

% Overall fit
%plot(x_values,std_fitted_all,'--','LineWidth',4,'color',[.5 .5 .5]);

% add legend
theLegendNames = unique(myPhosphoAuxiliary.myLegendNames, 'stable'); % b/c repetitions not in list of lines, get unique names.
figure(5); legend(lines4, theLegendNames,'Location','best');

h = findobj(gca,'Type','line')
x=get(h,'Xdata'),allx=x{:}
y=get(h,'Ydata'),ally=y{:}
xlim([0 max(allx)] )
ylim([0 max(ally)*1.05])

% Plot stds
% ===
fignr = 5;
fignr = fignr+1;
% Plot data

nrpositions = 3;
mymap = colorGray(nrpositions);
figure(fignr); clf; hold on;
set(gca,'FontSize',20);
%title('Growth speed vs. noise (1 point = 1 ''individual'')');
xlabel('\mu (doublings/hr)');
ylabel('std. dev.');
%axis([0 2 0 2])
lines4=[];

% Plot all 
% 732
lgh4 = plot(all_means_732,all_stds_732,'o','color',colorAmolfBlue)
lines4 = [lines4 lgh4];  % legend line
% 734
lgh4 = plot(all_means_734,all_stds_734,'o','color',colorAmolfYellow)
lines4 = [lines4 lgh4]; % legend line
% 735
lgh4 = plot(all_means_735,all_stds_735,'o','color',colorAmolfDarkGreen)
lines4 = [lines4 lgh4]; % legend line

%% plot length vs. speed
fignr = fignr+1;

figure(fignr); clf; hold on;
xlabel('length'); ylabel('growth speed (dbl/hr)');
title('WT speeds');
ylim([0,1.2]);

figure(fignr+1); clf; hold on;
xlabel('length'); ylabel('noise (dbl/hr)');
title('WT noises');
ylim([0,0.5]);

% WT rep1
posname='pos1crop';posdate='2014-05-01';
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir',myPhosphoAuxiliary.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none'),'quickMode',1);
figure(fignr); %plot([schnitzcells.length_fitNew],[schnitzcells.muP11_all],'x','color',colorAmolfBlue); % plot raw
xdata = [schnitzcells.length_fitNew]; ydata = [schnitzcells.muP11_all]; % load data
[the_bin_centers, the_means, the_stds, the_noises] = xy_mean_std_per_window(xdata, ydata, nr_windows, 1); % get stats
plot(the_bin_centers, the_means,'xk'); % plot stats in same graph
figure(fignr+1); plot(the_bin_centers, the_noises,'xk'); % plot noises in 2nd graph

% WT rep2
posname='pos4crop';posdate='2014_06_18';
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir',myPhosphoAuxiliary.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none'),'quickMode',1);
figure(fignr); %plot([schnitzcells.length_fitNew],[schnitzcells.muP11_all],'x','color',colorAmolfBlue); % plot raw
xdata = [schnitzcells.length_fitNew]; ydata = [schnitzcells.muP11_all]; % load data
[the_bin_centers, the_means, the_stds, the_noises] = xy_mean_std_per_window(xdata, ydata, nr_windows, 1); % get stats
plot(the_bin_centers, the_means,'ok'); % plot stats in same graph
figure(fignr+1); plot(the_bin_centers, the_noises,'ok'); % plot noises in 2nd graph

% WT rep3
posname='pos2crop';posdate='2014_06_18';
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir',myPhosphoAuxiliary.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none'),'quickMode',1);
figure(fignr); %plot([schnitzcells.length_fitNew],[schnitzcells.muP11_all],'x','color',colorAmolfBlue); % plot raw
xdata = [schnitzcells.length_fitNew]; ydata = [schnitzcells.muP11_all]; % load data
[the_bin_centers, the_means, the_stds, the_noises] = xy_mean_std_per_window(xdata, ydata, nr_windows, 1); % get stats

plot(the_bin_centers, the_means,'sk'); % plot stats in same graph
figure(fignr+1); plot(the_bin_centers, the_noises,'sk'); % plot noises in 2nd graph

fignr = fignr+2;

figure(fignr); clf; hold on;
xlabel('length'); ylabel('growth speed (dbl/hr)');
title('mutant speeds');
ylim([0,1.2]);

figure(fignr+1); clf; hold on;
xlabel('length'); ylabel('noise (dbl/hr)');
title('mutant noises');
ylim([0,0.5]);

% mut rep1
posname='pos8crop';posdate='2014-05-02';
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir',myPhosphoAuxiliary.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none'),'quickMode',1);
figure(fignr); %plot([schnitzcells.length_fitNew],[schnitzcells.muP11_all],'x','color',colorAmolfDarkGreen); % plot raw
xdata = [schnitzcells.length_fitNew]; ydata = [schnitzcells.muP11_all]; % load data
[the_bin_centers, the_means, the_stds, the_noises] = xy_mean_std_per_window(xdata, ydata, nr_windows, 1); % get stats
plot(the_bin_centers, the_means,'xk'); % plot stats in same graph
figure(fignr+1); plot(the_bin_centers, the_noises,'xk'); % plot noises in 2nd graph

% mut rep2
posname='pos2crop';posdate='2014_06_22';
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir',myPhosphoAuxiliary.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none'),'quickMode',1);
figure(fignr); %plot([schnitzcells.length_fitNew],[schnitzcells.muP11_all],'x','color',colorAmolfDarkGreen); % plot raw
xdata = [schnitzcells.length_fitNew]; ydata = [schnitzcells.muP11_all]; % load data
[the_bin_centers, the_means, the_stds, the_noises] = xy_mean_std_per_window(xdata, ydata, nr_windows, 1); % get stats
plot(the_bin_centers, the_means,'ok'); % plot stats in same graph
figure(fignr+1); plot(the_bin_centers, the_noises,'ok'); % plot noises in 2nd graph

% mut rep3
posname='pos3crop';posdate='2014_06_22';
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir',myPhosphoAuxiliary.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none'),'quickMode',1);
figure(fignr); %plot([schnitzcells.length_fitNew],[schnitzcells.muP11_all],'x','color',colorAmolfDarkGreen); % plot raw
xdata = [schnitzcells.length_fitNew]; ydata = [schnitzcells.muP11_all]; % load data
[the_bin_centers, the_means, the_stds, the_noises] = xy_mean_std_per_window(xdata, ydata, nr_windows, 1); % get stats
plot(the_bin_centers, the_means,'sk'); % plot stats in same graph
figure(fignr+1); plot(the_bin_centers, the_noises,'sk'); % plot noises in 2nd graph


%% divide it up in windows 
% (requires previous section)

xdata = [schnitzcells.length_fitNew]; ydata = [schnitzcells.muP11_all];
nr_windows = 50;

[the_bin_centers, the_means, the_stds, the_noises] = xy_mean_std_per_window(xdata, ydata, nr_windows);
%figure;
%plot(the_bin_centers, the_means,'xk');
%errorbar(the_bin_centers, the_means,the_noises,'xk');

figure;
plot(the_bin_centers, the_means,'xk');
figure;
plot(the_bin_centers, the_noises,'xk');






