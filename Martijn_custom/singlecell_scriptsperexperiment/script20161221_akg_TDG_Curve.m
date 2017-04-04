%%
%
% Plotting optimality curve from cAMP platereader dataset w. ASC1004,
% 2016_10_19 dataset.

%%
% Run first: 

if ~exist('DONTLOAD','var')
    error('Load dataset or set DONTLOAD..');
    % load('\\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\C_Platereader\2016_12_21_akg_tdg\23-Dec-2016CompleteAnalyzedData_GFP.mat');   
    % OR
    % load('\\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\C_Platereader\2017_01_03_akg_tdg_repeat2\05-Jan-2017CompleteAnalyzedData_GFP.mat');   
end
% See also \\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\C_Platereader\2016_12_21_akg_tdg
% OR \\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\C_Platereader\2017_01_03_akg_tdg_repeat2

%% Concentrations specific to this experiment
SHOWMANUALMUS = 0;

if ~exist('WHICHDATASET','var')
    WHICHDATASET = 1; % TDG
    %WHICHDATASET = 2; % AKG
end

if WHICHDATASET==1
    % prepare x-axis
    myprelConcentrations=[0 1.*.5.^[0 1 2 3 4]];
    myConcentrations=[myprelConcentrations];
elseif WHICHDATASET==2
    % prepare x-axis 
    myprelConcentrations=[0, 20 , 10, 5, 2.5];
    myConcentrations=[myprelConcentrations];
end

mean_ODplateaus = arrayfun(@(x) mean(output(x).ODPlateaus(:)),1:numel(output));
std_ODplateaus = arrayfun(@(x) std(output(x).ODPlateaus(:)),1:numel(output));

plateauData = [mean_ODplateaus];

%% Acquire mu values
% Run ExtractFitPlateReaderData_General_Part3_Plotting first

muData = [output(:).muValuesMean];

%% Create dual y-axis optimality plot, 

hF = figure; clf; 

% plot individual well points
% ....
% manualMuValues

% plot lines
indicesNotZero = ~myConcentrations==0;
[ax,l1,l2]  = plotyy(myConcentrations(indicesNotZero),muData(indicesNotZero),...
                     myConcentrations(indicesNotZero),plateauData(indicesNotZero));%,'semilogx');
set(l1,'Marker','o','LineWidth',3,'Color','r');
set(l2,'Marker','s','LineWidth',2,'Color',[.7 .7 .7]);


hold(ax(1), 'on');
hold(ax(2), 'on');

% x lim
myXlim = [min(myConcentrations(myConcentrations>0))/3,max(myConcentrations)*3];
xlim(ax(1),myXlim);
xlim(ax(2),myXlim);

% plot values at 0
l0 = plot(ax(1),myXlim(1),muData(~indicesNotZero));
set(l0,'Marker','o','LineWidth',3,'Color','r');
l0 = plot(ax(2),myXlim(1),plateauData(~indicesNotZero));
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
ylim(ax(2),[0,max(plateauData)*2]);

MW_makeplotlookbetter(15);

if WHICHDATASET==1
    xlabel('Concentration TDG [mM]','Color','k');
elseif WHICHDATASET ==2
    xlabel('Concentration aKG [mM]','Color','k');
end


ylabel(ax(1),'Growth rate [dbl/hr]','Color','k');
ylabel(ax(2),'OD plateau value [a.u.]','Color','k');
set(ax(1), 'YColor', 'k');
set(ax(2), 'YColor', 'k');

legend([l1,l2],{'Growth rates','Max. OD value observed'})

% Per default larger fonts don't fit this window
set(ax(1), 'Position',[0.15 0.15 0.65 0.8]);

% Set ticks
labelLocations=logspace(1,4,4);
%set(gca,'XTick',labelLocations);%'XTickLabel',USERSETTINGS.wellNamesToPlot, 
set(ax(1),'YTick',[0:0.2:1]);%'XTickLabel',USERSETTINGS.wellNamesToPlot, 
set(ax(2),'YTick',[0:0.1:max(plateauData)*2]);%'XTickLabel',USERSETTINGS.wellNamesToPlot, 



%%

%OUTPUTDIR = '\\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\C_Platereader\2016_12_21_akg_tdg\'

if ~exist('OUTPUTDIR','var')
    disp('set OUTPUTDIR to output plots');
else
    names={'TDG','aKG'};
    saveas(hF, [OUTPUTDIR 'curve_' names{WHICHDATASET} '.fig']);
    saveas(hF, [OUTPUTDIR 'curve_' names{WHICHDATASET} '.svg']);
    saveas(hF, [OUTPUTDIR 'curve_' names{WHICHDATASET} '.tif']);
end




