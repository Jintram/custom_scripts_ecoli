
%% Settings

some_colors;

% Common settings
myRootDir = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\'

%% Load data_______________________________________________________________

% 735___
p = DJK_initschnitz('pos8crop','2014-05-02','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

[plotx_735_pos8, ploty_735_pos8] = ...
    DJK_plot_scatterColor(p, s_all, 'av_mu_fitNew', 'av_time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);

% 732___
p = DJK_initschnitz('pos1crop','2014-05-01','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

[plotx_732_pos1, ploty_732_pos1] = ...
    DJK_plot_scatterColor(p, s_all, 'av_mu_fitNew', 'av_time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);


%% Plotting the data________________________________________________________

% Scatter plot
% ====

figure(1);
clf(1);
hold on;


h732_pos1 = scatter(plotx_732_pos1, ploty_732_pos1, 100,colorAmolfGreen,'o','LineWidth',2); % WT, 37C
h735_pos8 = scatter(plotx_735_pos8, ploty_735_pos8, 100,colorAmolfBlue,'o','LineWidth',2); % WT, 37C

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







