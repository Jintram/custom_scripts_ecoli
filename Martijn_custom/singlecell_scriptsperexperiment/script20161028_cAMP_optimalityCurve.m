%%
%
% Plotting optimality curve from cAMP platereader dataset w. ASC1004,
% 2016_10_19 dataset.

PLOTMANUALFITS=1;

%%
% Run first: 
if ~exist('alreadyLoaded','var')
    %load('U:\ZZ_EXPERIMENTAL_DATA\Platereader\2016_10_19_asc1004_cAMPseries_part1\28-Oct-2016CompleteAnalyzedData_GFP.mat');
    load 'U:\ZZ_EXPERIMENTAL_DATA\Platereader\2016_10_19_asc1004_cAMPseries_part1\22-Mar-2017CompleteAnalyzedData_GFP.mat';
    alreadyLoaded=1; 
    USERSETTINGS.plotManualFits=PLOTMANUALFITS; % note this parameter is first loaded in load prev. line
    
    USERSETTINGS.wellNamesToPlot= {'asc1004A','asc1004B', 'asc1004C', 'asc1004D', 'asc1004E', 'asc1004F', 'asc1004G', 'asc1004H', 'asc1004zero'}
    ExtractFitPlateReaderData_General_Part3_Plotting        
end

% cosmetics
some_colors;

%{
%OR REDO ANALYSIS:
USERSETTINGS.myRootDir='U:\ZZ_EXPERIMENTAL_DATA\Platereader\';
USERSETTINGS.myDateDir='2016_10_19_asc1004_cAMPseries_part1\';
USERSETTINGS.datafile= '2016_10_19_cAMPDilutionASC1004_part1.xls';
USERSETTINGS.customSuffix = '_OD';
USERSETTINGS.ODorFluor = 1;

USERSETTINGS.ODmin=0.05;
USERSETTINGS.ODmax=0.10;

USERSETTINGS.fitManual = 1;

% TIMEINDEXES=[7,9,11], YINDEXES   = [8,10,12] % new platereader
TIMEINDEXES=[5,7,9], YINDEXES   = [6,8, 10] % old platereader

ExtractFitPlateReaderData_General_Part1
ExtractFitPlateReaderData_General_Part2_OD

% GFP
USERSETTINGS.myRootDir='U:\ZZ_EXPERIMENTAL_DATA\Platereader\';
USERSETTINGS.myDateDir='2016_10_19_asc1004_cAMPseries_part1\';
USERSETTINGS.datafile= '2016_10_19_cAMPDilutionASC1004_part1.xls';
USERSETTINGS.customSuffix = '_GFP';
USERSETTINGS.ODorFluor = 2;
USERSETTINGS.platereader = 'OLD';

% for fluor: only if you want to REDO the manual range based on fluor
USERSETTINGS.fitManual = 0; 

%TIMEINDEXES=[5], YINDEXES   = [6] % new platereader
TIMEINDEXES=[11], YINDEXES   = [12] % old platereader

ExtractFitPlateReaderData_General_Part1
ExtractFitPlateReaderData_General_Part2_Fluor
 

USERSETTINGS.plotManualFits=1;

USERSETTINGS.wellNamesToPlot= {'asc1004A','asc1004B', 'asc1004C', 'asc1004D', 'asc1004E', 'asc1004F', 'asc1004G', 'asc1004H', 'asc1004zero'}

ExtractFitPlateReaderData_General_Part3_Plotting
%}

%% Concentrations specific to this experiment

myprelConcentrations=[10000./3.33.^[0,1,2,3,4,5,6,7]];
myConcentrations=[myprelConcentrations(1:7)+(1/601*800) myprelConcentrations(8)+(1/857.5*800)];

%% Create dual y-axis optimality plot, 
if ~isfield(USERSETTINGS,'plotManualFits')
        USERSETTINGS.plotManualFits = 1;
end

optimH = figure(1); clf; 

% plot individual well points
% ....
% manualMuValues

% plot lines
if USERSETTINGS.plotManualFits
    muData = [output.manualMuValuesMean];
else
    muData = [output.muValuesMean];
end
plateauData = arrayfun(@(x) mean(output(x).ODPlateaus), 1:numel(output)); %[output.ODPlateaus];
[ax,l1,l2]  = plotyy(myConcentrations,muData(1:8),...
                     myConcentrations,plateauData(1:8));%,'semilogx');
set(l1,'Marker','o','LineWidth',3,'Color',colorblind(2,:));
set(l2,'Marker','s','LineWidth',2,'Color',[.7 .7 .7]);


hold(ax(1), 'on');
hold(ax(2), 'on');

% x lim
myXlim = [min(myConcentrations(myConcentrations>0))/3,max(myConcentrations)*3];
xlim(ax(1),myXlim);
xlim(ax(2),myXlim);

% plot values at 0
l0 = plot(ax(1),myXlim(1),muData(end));
set(l0,'Marker','o','LineWidth',3,'Color',colorblind(2,:));
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

myYMax=max(muData)*1.4;
ylim(ax(1),[0,myYMax]);
ylim(ax(2),[min(plateauData),max(plateauData)*1.4]);

MW_makeplotlookbetter(15);

xlabel('Concentration cAMP [uM]','Color','k');
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

%% Plot w/o max OD, 2 plots for manual & range
optimH2 = figure(2); clf; hold on;

muDataManual = [output.manualMuValuesMean];
muDataRange = [output.muValuesMean];

l1=plot(myConcentrations,muDataManual(1:8));
set(l1,'Marker','o','LineWidth',3,'Color',colorblind(2,:));
l2=plot(myConcentrations,muDataRange(1:8));
set(l2,'Marker','o','LineWidth',3,'Color',colorblind(3,:));

ylabel('Growth rate [dbl/hr]','Color','k');
xlabel('Concentration cAMP [uM]','Color','k');

hold on;

% x lim
myXlim = [min(myConcentrations(myConcentrations>0))/3,max(myConcentrations)*3];
xlim(myXlim);

% plot values at 0
l0 = plot(myXlim(1),muDataManual(end));
set(l0,'Marker','o','LineWidth',3,'Color',colorblind(2,:));
l0 = plot(myXlim(1),muDataRange(end));
set(l0,'Marker','o','LineWidth',3,'Color',colorblind(3,:));

% log scale x axes
set(gca,'xscale','log');

myYMax=max(muData)*1.4;
ylim(ax(1),[0,myYMax]);

% Per default larger fonts don't fit this window
%set(ax(1), 'Position',[0.15 0.15 0.65 0.8]);

% Set ticks
labelLocations=logspace(1,4,4);
set(gca,'XTick',labelLocations);%'XTickLabel',USERSETTINGS.wellNamesToPlot, 
set(gca,'YTick',[0:0.2:1]);%'XTickLabel',USERSETTINGS.wellNamesToPlot, 

MW_makeplotlookbetter(20);

%% Plot Benjamin's data in here
figure(optimH2);

concentrationsBT = [10000./2.^[0 1 2 3 4 5 6 7 8 9]];
dataBT =        [NaN    0.2662    0.4116    0.5112    0.5353    0.5513    0.5289    0.3775    0.2384       NaN; ...
                 NaN    0.2541    0.4316    0.5158    0.5384    0.5601    0.5242    0.3547    0.2072       NaN; ...
                 0.0659    0.2694    0.4541    0.5198    0.5549    0.5587    0.5218    0.3573    0.2294       NaN];
dataBTbase2 = dataBT/log(2);
             
meandataBT = mean(dataBT);
meandataBTbase2 = mean(dataBTbase2);
stddataBT = std(dataBT);

l3=plot(concentrationsBT,meandataBTbase2);
set(l3,'Marker','o','LineWidth',3,'Color',colorblind(4,:));

%legend([l1,l2,l3],{'Growth rates','Max. OD value observed','Data Benjamin'})
%legend([l1,l3],{'Data Martijn','Data Benjamin'})
legend([l1,l2,l3],{'Data Martijn (method 1)','Data Martijn (method 2)','Data Benjamin'})

myYMax=max([muData meandataBTbase2])*1.4;
ylim([0,myYMax]);

MW_makeplotlookbetter(20);

%% Plot current low/opt/high concentrations in here

myExternalCampConcentrations=[80,800,5000];

for i=1:numel(myExternalCampConcentrations)
    %plot([myExternalCampConcentrations(i),myExternalCampConcentrations(i)],[0,myYMax],'--k','LineWidth',2)
    plot([myExternalCampConcentrations(i)],[0],'^k','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10)
end

%% Make a plot of BT data for exporting the line

figure(3); clf; 

subplot(1,2,1); hold on;
plot(concentrationsBT,meandataBTbase2,'-ok','LineWidth',3);

subplot(1,2,2); hold on;
plot(concentrationsBT,meandataBTbase2,'-ok','LineWidth',3);
set(gca,'xscale','log');



