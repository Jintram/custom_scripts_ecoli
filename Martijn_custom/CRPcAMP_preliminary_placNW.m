% CRP.cAMP project
% ===
% MW 2015/03/27

% Plot some of Noreen's extr/intr noise data, but with different aim.
% Data is available at:
% F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal
% Data is extracted from:
% \\biofysicasrv\Users2\Walker\ExtrinsicNoise\Data_Collection\Data
% See for more info:
% \\biofysicasrv\Users2\Walker\ExtrinsicNoise\Data_Collection\List_Experiments_Used_and_Unsused_ExtrNoiseProject.xlsx
% (Note that 1/3 of cellcycle means N datapoints used, which can be seen
% from name mu parameter.)

%% Load dataset
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Maltose basal 2012-05-16 pos1.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Maltose basal 2012-07-26 pos4.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Acetate basal 2012-07-19 pos4.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\RDM basal 2011-10-19 pos1.mat');
load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\RDM basal 2011-10-19 pos3.mat');

% clear all, close all

%% Simple params
cyan    =  [0,1,1]/2;
yellow  =  [212/255, 170/255, 0];

growthRateRange = [prctile(mu_Third_cycCor,5), prctile(mu_Third_cycCor,95)];
fullGrowthRateRange = [min(mu_Third_cycCor), max(mu_Third_cycCor)];

%% Plot histogram

% Plot 
figure, hold on
[histydC,histxdC] = hist(dC5_cycCor,100)
[histydY,histxdY] = hist(dY5_cycCor,100)
histxdC=histxdC/sum(histxdC);
histxdY=histxdY/sum(histxdY);
plot(histxdC,histydC,'-','Color',cyan, 'Linewidth',3), hold on
plot(histxdY,histydY,'-','Color',yellow, 'Linewidth',3), hold on
xlabel('Production rate')
ylabel('Normalized count')
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal')
set(gca,'FontSize',15)


%% Plot scatter of production rates
% settings
plotpretty = 0

% plotting
hFig = figure
offset=100; width=500; height=500;
set(hFig, 'Position', [offset offset width height])

if plotpretty
    scatter_patches(mu_Third_cycCor,dC5_cycCor,10,'o','FaceColor',cyan,'EdgeColor','none','FaceAlpha',0.3);
    scatter_patches(mu_Third_cycCor,dY5_cycCor,10,'o','FaceColor',yellow,'EdgeColor','none','FaceAlpha',0.3);
else
    subplot(2,1,1), hold on
    plot(mu_Third_cycCor,dC5_cycCor,'o','Color',cyan);
    subplot(2,1,2), hold on
    plot(mu_Third_cycCor,dY5_cycCor,'o','Color',yellow);
end

for i = 1:2
    subplot(2,1,i)
    xlabel('Growth rate (dbl/hr)')
    ylabel('Producation rate')
    %Set all fontsizes
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal')
    set(gca,'FontSize',15)
    xlim(fullGrowthRateRange)
    %ylim([prctile(dC5_cycCor,5), prctile(dC5_cycCor,95)])
    disp('Note some data falls outside range')
end

% fit data and plot
pdC = polyfit(mu_Third_cycCor',dC5_cycCor',1);
fitydC = pdC(1)*growthRateRange + pdC(2);
pdY = polyfit(mu_Third_cycCor',dY5_cycCor',1);
fitydY = pdY(1)*growthRateRange + pdY(2);

subplot(2,1,1)
plot(growthRateRange,fitydC,'-','Color','k','LineWidth',3)
subplot(2,1,2)
plot(growthRateRange,fitydY,'-','Color','k','LineWidth',3)



%mydC5_cycCor-dC5_cycCor

%% Extract the concentration data using Daan's function

% Prepare dummy p struct
p={}, p.movieName = 'dummy';
% Retrieve data C5_mean_cycCor
retrievedData = DJK_get_schnitzData(p, schnitzcells_rm, 'time_atY','dataFields',{'C5_mean_cycCor'},'fitTime',fitTime);
C5_mean_cycCor = [retrievedData.C5_mean_cycCor];
% Retrieve data Y5_mean_cycCor
retrievedData = DJK_get_schnitzData(p, schnitzcells_rm, 'time_atY','dataFields',{'Y5_mean_cycCor'},'fitTime',fitTime);
Y5_mean_cycCor = [retrievedData.Y5_mean_cycCor];



%% Plot scatter of production rates
% Settings
prettyplotting=0


% Plotting
hFig = figure
offset=100; width=500; height=500;
set(hFig, 'Position', [offset offset width height])

if prettyplotting
    scatter_patches(mu_Third_cycCor,C5_mean_cycCor,10,'o','FaceColor',cyan,'EdgeColor','none','FaceAlpha',0.3);
    scatter_patches(mu_Third_cycCor,Y5_mean_cycCor,10,'o','FaceColor',yellow,'EdgeColor','none','FaceAlpha',0.3);
else
    subplot(2,1,1), hold on
    plot(mu_Third_cycCor,C5_mean_cycCor,'x','Color',cyan);
    subplot(2,1,2), hold on
    plot(mu_Third_cycCor,Y5_mean_cycCor,'x','Color',yellow);
end

for i = 1:2
    subplot(2,1,i)
    xlabel('Growth rate (dbl/hr)')
    ylabel('Concentration')
    %Set all fontsizes
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal')
    set(gca,'FontSize',15)

    xlim(fullGrowthRateRange)
    %ylim([-200,800])
    disp('Note some data falls outside range')
end

% fit data and plot
pC = polyfit(mu_Third_cycCor',C5_mean_cycCor',1);
fityC = pC(1)*growthRateRange + pC(2);
pY = polyfit(mu_Third_cycCor',Y5_mean_cycCor',1);
fityY = pY(1)*growthRateRange + pY(2);

subplot(2,1,1)
plot(growthRateRange,fityC,'-','Color','k','LineWidth',3)
subplot(2,1,2)
plot(growthRateRange,fityY,'-','Color','k','LineWidth',3)





%% Sanity check

% Check if I retrieve same data as Noreen.
retrievedData = DJK_get_schnitzData(p, schnitzcells_rm, 'time_atY','dataFields',{'dC5_cycCor'},'fitTime',fitTime);
mydC5_cycCor = [retrievedData.dC5_cycCor];
sanityCheck = mydC5_cycCor - dC5_cycCor;
sanityCheck(find(sanityCheck>0))
disp('sanityCheck chould be empty');

%{ 
Old code

%% Collect data from schnitzcells_rm manually
% It's more convenient to use DJK_get_schnitzData.
% Collect data
% Select instead of: myYdata = [schnitzcells_rm(:).Y5_mean]
mydC5_cycCor = [];
for i = 1:numel(schnitzcells_rm)
    if schnitzcells_rm(i).useForPlot == 1 & ...
            schnitzcells_rm(i).time > fitTime(1) & ...
            schnitzcells_rm(i).time < fitTime(2)
        mydC5_cycCor = [mydC5_cycCor schnitzcells_rm(i).dC5_cycCor];
    end
end
mydC5_cycCor

%}



