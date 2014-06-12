














% Newer version available, look at phospho_base_noise.m now!
























%% Settings

some_colors;

% Common settings
myRootDir = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\'

%% Load data_______________________________________________________________

% 735___
p = DJK_initschnitz('pos8crop','2014-05-02','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

% select which Schnitzcells to take into account (all)
s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

[plotx_735_pos8, ploty_735_pos8] = ...
    DJK_plot_scatterColor(p, s_all, 'muP11_all', 'time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 1);

% 733___
p = DJK_initschnitz('pos8crop','2014-05-01','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

% select which Schnitzcells to take into account (all)
s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

[plotx_733_pos8, ploty_733_pos8] = ...
    DJK_plot_scatterColor(p, s_all, 'muP11_all', 'time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 1);


% 732___
p = DJK_initschnitz('pos1crop','2014-05-01','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

[plotx_732_pos1, ploty_732_pos1] = ...
    DJK_plot_scatterColor(p, s_all, 'muP11_all', 'time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 1);

   
%% Plotting the data________________________________________________________

% Note that this is also performed by Daan's script. But this code puts
% the data in the same plot.

% Scatter plot
% ====

figure(1);
clf(1);
hold on;


h732_pos1 = scatter(plotx_732_pos1, ploty_732_pos1, 100,colorAmolfGreen,'x','LineWidth',1); % WT, 37C
h735_pos8 = scatter(plotx_735_pos8, ploty_735_pos8, 100,colorAmolfBlue,'x','LineWidth',1); % WT, 37C

%h42_pos6 = scatter(plotx_37_pos6, ploty_37_pos6, 100,colorAmolfBlue,'x','LineWidth',2); % 717, 37C
%h2 = scatter(plotx_pos2, ploty_pos2, 100,'r','o','LineWidth',2); % WT, 42C
%h13 = scatter(plotx_pos13, ploty_pos13, 100,'r','x','LineWidth',2); % 717, 42C

legend('732','735');
%legend('WT, 37C','717, 37C','WT, 42C','717, 42C');

%set(h13)

set(gca,'FontSize',20);
title('Doubling time for each ''individual''');
xlabel('time (min)');
ylabel('\mu (doublings/hr)');

hold off;

    
%% Preparing histograms
% =====

% bins
dx=2/50;
mybins=[0:dx:2];

% hist
[score732_pos1] = hist(ploty_732_pos1,mybins);
[score735_pos8] = hist(ploty_735_pos8,mybins);

% normalizing
A732=sum(score732_pos1)*dx;
A735=sum(score735_pos8)*dx;

score732_pos1=score732_pos1./A732;
score735_pos8=score735_pos8./A735;

%% plotting
figure(2);
clf(2);
hold on;

plot(mybins, score732_pos1,'--o','LineWidth',3,'color',colorAmolfGreen); 
plot(mybins, score735_pos8,'--o','LineWidth',3,'color',colorAmolfBlue); 

legend('732','735');

set(gca,'FontSize',20);
title('Distribution of doubling time');
xlabel('\mu (doublings/hr)');
ylabel('Probability (normalized)');

hold off;


%% Sliding window variance

nrpositions = 3;
mylegendnames={};
mymap = colorGray(nrpositions);
fignumstart=1;
figure(fignumstart); clf; figure(fignumstart+1); clf; figure(fignumstart+2); clf;

% WT
mylegendnames = [mylegendnames, 'Wildtype', 'WT fit'];
datax=plotx_732_pos1;
datay=ploty_732_pos1;
sliding_window_stds_plot(datax,datay,fignumstart,mymap(1,:),1,1);

% 733
mylegendnames = [mylegendnames, '\Delta{pykF} mutant ', 'mutant fit'];
datax=plotx_733_pos8;
datay=ploty_733_pos8;
sliding_window_stds_plot(datax,datay,fignumstart,mymap(2,:),1,1);

% 735
mylegendnames = [mylegendnames, '\Delta{pykF}/\Delta{ppc} mutant ', 'mutant fit'];
datax=plotx_735_pos8;
datay=ploty_735_pos8;
sliding_window_stds_plot(datax,datay,fignumstart,mymap(2,:),1,1);

figure(fignumstart); legend(mylegendnames(1:2:length(mylegendnames)),'Location','best');
figure(fignumstart+1); legend(mylegendnames,'Location','best');
figure(fignumstart+2); legend(mylegendnames,'Location','best');



