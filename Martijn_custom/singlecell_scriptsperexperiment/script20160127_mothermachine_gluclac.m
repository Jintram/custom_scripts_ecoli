





OUTPUTDIR = 'D:\Local_Software\Martijn_extensions\Martijn_custom\singlecell_scriptsperexperiment\plots20160127_mothermachine_gluclac\';
FILEPREFIX = [p.movieDate p.movieName '_'];



FIELDSTOPLOT      = {'Y6_mean_all',     'C6_mean_all'};     %,      'dY5_cycCor',       'dC5_cycCor'};
TIMEFIELDS        = {'time_atY',        'time_atC'};        %,         'time_atdY',         'time_atdC'};
CORRESPONDINGLIMS = {[0 600],           [0 200]};           %,            [],                 []};
YLABELNAMES       = {'CRP reporter',    'Const. reporter'}; %,  'Production CRP',   'Production const. reporter'};

% Time data
% Part 1
% 13-01-2017 19:00 Lac
% 13-01-2017 24:00 Gluc
% 14-01-2017 05:00 Lac
% 14-01-2017 10:00 Gluc
% 14-01-2017 15:00 Lac
% 14-01-2017 20:00 Gluc
% Part 2 (coincidentally exactly following previous pattern)
% 15-01-2017 01:00 Lac
% 15-01-2017 06:00 Gluc
% 15-01-2017 11:00 Lac
% 15-01-2017 16:00 Gluc
% 15-01-2017 21:00 Lac
% 15-01-2017 02:00 Gluc

dataY={};
for plotIndex = 1:numel(FIELDSTOPLOT)

    if ~isempty(CORRESPONDINGLIMS{plotIndex})
        bottomlim   =CORRESPONDINGLIMS{plotIndex}(1);
        toplim      =CORRESPONDINGLIMS{plotIndex}(2);
    else
        bottomlim=0;
        toplim=0;
    end

    blockSignalTime    = [-5         0           0       5       5          10        10     15     15        20        20      25];
    blockSignalHighLow = [bottomlim  bottomlim   toplim  toplim  bottomlim  bottomlim toplim toplim bottomlim bottomlim toplim  toplim];

    dataX = [schnitzcells.(TIMEFIELDS{plotIndex})];
    dataXHour = dataX./60;
    %dataY = schnitzcells.time_atC
    dataY{plotIndex} = [schnitzcells.(FIELDSTOPLOT{plotIndex})];

    h=figure; clf; hold on;

    plot(blockSignalTime, blockSignalHighLow,'-','Color',[.5 .5 .5],'LineWidth',3);

    %plot(dataXHour, dataY,'.');
    scatter(dataXHour, dataY{plotIndex},7^2,'filled', 'MarkerFaceAlpha',.1);

    % Calculate average points
    timePoints=unique(dataXHour);
    dT = timePoints(2)-timePoints(1);
    binsX = [timePoints timePoints(end)+dT] - dT/2;
    [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=binnedaveraging({dataXHour},{dataY{plotIndex}},binsX);
    plot(binCenters,meanValuesForBins,'ok','MarkerFaceColor','k','MarkerSize',10);
   
    if ~isempty(CORRESPONDINGLIMS{plotIndex})
        ylim(CORRESPONDINGLIMS{plotIndex});
    end
    xlim([-5,25]);

    xlabel('time [hrs]');
    ylabel(YLABELNAMES{plotIndex});

    MW_makeplotlookbetter(20);

    saveas(h, [OUTPUTDIR 'SVG_' FILEPREFIX FIELDSTOPLOT{plotIndex} '.svg']);
    saveas(h, [OUTPUTDIR 'TIF_' FILEPREFIX FIELDSTOPLOT{plotIndex} '.tif']);
    saveas(h, [OUTPUTDIR 'FIG_' FILEPREFIX FIELDSTOPLOT{plotIndex} '.fig']);
    
end

%%

h=figure; clf; hold on;

xlabel('time [hrs]');
ylabel('CRP label/const. label');

MW_makeplotlookbetter(20);

% blocksignal
toplim=CORRESPONDINGLIMS{1}(2)/CORRESPONDINGLIMS{2}(2)*2;
bottomlim=0;
blockSignalTime    = [-5         0           0       5       5          10        10     15     15        20        20      25];
blockSignalHighLow = [bottomlim  bottomlim   toplim  toplim  bottomlim  bottomlim toplim toplim bottomlim bottomlim toplim  toplim];
plot(blockSignalTime, blockSignalHighLow,'-','Color',[.5 .5 .5],'LineWidth',3);

% data
normalizedData = dataY{1}./dataY{2};
scatter(dataXHour, normalizedData,7^2,'filled', 'MarkerFaceAlpha',.1);

% Calculate average points
timePoints=unique(dataXHour);
dT = timePoints(2)-timePoints(1);
binsX = [timePoints timePoints(end)+dT] - dT/2;
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=binnedaveraging({dataXHour},{normalizedData},binsX);
plot(binCenters,meanValuesForBins,'ok','MarkerFaceColor','k','MarkerSize',10);
plot(binCenters,meanValuesForBins,'-k','LineWidth',2);

ylim([0, toplim]);
xlim([-5,25]);

saveas(h, [OUTPUTDIR 'SVG_' FILEPREFIX 'CRPdivConst.svg']);
saveas(h, [OUTPUTDIR 'TIF_' FILEPREFIX 'CRPdivConst.tif']);
saveas(h, [OUTPUTDIR 'FIG_' FILEPREFIX 'CRPdivConst.fig']);



