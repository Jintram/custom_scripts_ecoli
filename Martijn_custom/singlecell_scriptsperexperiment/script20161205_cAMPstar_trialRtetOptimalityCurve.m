%%
%
% Plotting optimality curve from cAMP platereader dataset w. ASC1004,
% 2016_10_19 dataset.

%%
% Run first: 

if ~exist('DONTLOAD','var')
    error('Load dataset or set DONTLOAD..');
    % load('U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\2016_12_05_CRPstar_Rtet\05-Dec-2016CompleteAnalyzedData_GFP.mat');
end
% See also U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\2016_12_05_CRPstar_Rtet

%% Concentrations specific to this experiment
SHOWMANUALMUS = 1;

% WHICHDATASET = 2;

if WHICHDATASET==1
    % prepare x-axis
    myprelConcentrations=[100./3.33.^[0,1,2,3,4,5]];
    myConcentrations=[myprelConcentrations 0];
elseif WHICHDATASET==2
    % prepare x-axis 
    myprelConcentrations=[200./3.33.^[0,1,2,3,4,5 6]];
    myConcentrations=[myprelConcentrations 0];
end

% prepare data (since we have no 0 this time, I add a NaN value
if SHOWMANUALMUS
    muData = [output.manualMuValuesMean NaN];
else
    muData = [output.muValuesMean NaN];
end

mean_ODplateaus = arrayfun(@(x) mean(output(x).ODPlateaus(:)),1:numel(output));
std_ODplateaus = arrayfun(@(x) std(output(x).ODPlateaus(:)),1:numel(output));

plateauData = [mean_ODplateaus NaN];

%% Create dual y-axis optimality plot, 

figure; clf; 

% plot individual well points
% ....
% manualMuValues

% plot lines
[ax,l1,l2]  = plotyy(myConcentrations(1:end-1),muData(1:end-1),...
                     myConcentrations(1:end-1),plateauData(1:end-1));%,'semilogx');
set(l1,'Marker','o','LineWidth',3,'Color','r');
set(l2,'Marker','s','LineWidth',2,'Color',[.7 .7 .7]);


hold(ax(1), 'on');
hold(ax(2), 'on');

% x lim
myXlim = [min(myConcentrations(myConcentrations>0))/3,max(myConcentrations)*3];
xlim(ax(1),myXlim);
xlim(ax(2),myXlim);

% plot values at 0
l0 = plot(ax(1),myXlim(1),muData(end));
set(l0,'Marker','o','LineWidth',3,'Color','r');
l0 = plot(ax(2),myXlim(1),plateauData(end));
set(l0,'Marker','s','LineWidth',2,'Color',[.7 .7 .7]);

% Put graph 1 on the foreground
% (Thanks to https://nl.mathworks.com/matlabcentral/answers/21181-setting-order-and-transparency-in-plotyy)
uistack(ax(1));
set(ax(1), 'Color', 'none');
set(ax(2), 'Color', 'w');
          
% log scale x axes
set(ax(1),'xscale','log');
set(ax(2),'xscale','log');

% Set cosmetics
set(ax(1),'fontsize',15);
set(ax(2),'fontsize',15);

ylim(ax(1),[0,1]);
ylim(ax(2),[0,0.5]);

MW_makeplotlookbetter(15);

xlabel('Concentration AtC [ng/mL]','Color','k');
ylabel(ax(1),'Growth rate [dbl/hr]','Color','k');
ylabel(ax(2),'OD plateau value [a.u.]','Color','k');
set(ax(1), 'YColor', 'k');
set(ax(2), 'YColor', 'k');

legend([l1,l2],{'Growth rates','Max. OD value observed'})

% Per default larger fonts don't fit this window
set(ax(1), 'Position',[0.15 0.15 0.65 0.8]);

% Set ticks
labelLocations=logspace(1,4,4);
set(gca,'XTick',labelLocations);%'XTickLabel',USERSETTINGS.wellNamesToPlot, 
set(ax(1),'YTick',[0:0.2:1]);%'XTickLabel',USERSETTINGS.wellNamesToPlot, 
set(ax(2),'YTick',[0:0.1:0.5]);%'XTickLabel',USERSETTINGS.wellNamesToPlot, 










