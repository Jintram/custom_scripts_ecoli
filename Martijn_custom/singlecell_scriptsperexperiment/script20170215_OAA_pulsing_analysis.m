
%% 

% 

% ACS990 (WT) data:
%{
load('H:\EXPERIMENTAL_DATA_2017\2017-02-02_OAA_pulsing_3\pos2crop\data\pos2crop-Schnitz.mat');
FRAMEDELAY=18/60;
TITLEPART2 = ', wild type';
run('H:\EXPERIMENTAL_DATA_2017\2017-02-02_OAA_pulsing_3\blockpulsedescriptionASC990_20170202.m')
%}


% ASC1004 (delta cAMP) data:
%{
% Note that this dataset has two subpopulations!
load('H:\EXPERIMENTAL_DATA_2017\2017-02-10_OAA_pulsing_asc1004\pos1crop\data\pos1crop-Schnitz.mat');
FRAMEDELAY = 319/60;
TITLEPART2 = ', \deltacAMP';
run('H:\EXPERIMENTAL_DATA_2017\2017-02-10_OAA_pulsing_asc1004\blockpulsedescription.m')
%}
%{
load('H:\EXPERIMENTAL_DATA_2017\2017-02-10_OAA_pulsing_asc1004\pos2crop\data\pos2crop-Schnitz.mat');
FRAMEDELAY = 1038/60; % For 2nd dataset
TITLEPART2 = ', \DeltacAMP';
run('H:\EXPERIMENTAL_DATA_2017\2017-02-10_OAA_pulsing_asc1004\blockpulsedescription.m')
%}





%% 

FIELDSTOTTAKE = {'muP5_fitNew'};
TIMEFIELDS ={'time_atC'};
TITLEPART1 = {'Growth dynamics'};

%{
FIELDSTOTTAKE = {'dY5_cycCor'};
TIMEFIELDS ={'time_atdY'};
TITLEPART1 = {'Fluor prod. rate'};
%}


%{
FIELDSTOTTAKE = {'Y5_mean', 'C5_mean'};
TITLEPART1 = {'CRP dynamics','Const. dynamics'};
%}
%%
for fieldIndex = 1:numel(FIELDSTOTTAKE)
    %% 
    currentFieldName = FIELDSTOTTAKE{fieldIndex};
    currentTimeField = TIMEFIELDS{fieldIndex};
    
    %%
    timeData = [schnitzcells.(currentTimeField)]/60;
    fluorData = [schnitzcells.(currentFieldName)];
    times=unique(timeData);
    if ~exist('MYYLIM','var')
        MYYLIM=[0,max(fluorData)*1.05];
    end

    figure(fieldIndex); clf; hold on;

    % All data
    scatter(timeData, fluorData,3^2,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor',[.5 .5 .5]);

    % Get and plot means
    for tprobe = times
        indices = timeData==tprobe;
        selectedTime = timeData(indices);
        selectedFluor = fluorData(indices);
        plot(mean(selectedTime),mean(selectedFluor),'ok','MarkerFaceColor','k')
    end

    % Plot block pulse
    ValuesToPlotRescaled=ValuesToPlot;
    ValuesToPlotRescaled(ValuesToPlot==1)=MYYLIM(2);
    ValuesToPlotRescaled(ValuesToPlot==2)=MYYLIM(1);%*.9542;
    plot(TimesToPlot-FRAMEDELAY,ValuesToPlotRescaled,'LineWidth',3,'Color','b');

    % Cosmetics
    xlim([0,42]);
    ylim(MYYLIM);

    title([TITLEPART1{fieldIndex} TITLEPART2]);
    MW_makeplotlookbetter(20);
end

%%
%{
timeData = [schnitzcells.time]/60;
fluorData = [schnitzcells.C6_mean];
times=unique(timeData);
%if ~exist('MYYLIM','var')
MYYLIM=[0,max(fluorData)*1.05];
%end

figure(2); clf; hold on;

% All data
scatter(timeData, fluorData,3^2,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor',[.5 .5 .5]);

% Get and plot means
for tprobe = times
    indices = timeData==tprobe;
    selectedTime = timeData(indices);
    selectedFluor = fluorData(indices);
    plot(mean(selectedTime),mean(selectedFluor),'ok','MarkerFaceColor','k')
end

% Plot block pulse
ValuesToPlotRescaled=ValuesToPlot;
ValuesToPlotRescaled(ValuesToPlot==1)=MYYLIM(2);
ValuesToPlotRescaled(ValuesToPlot==2)=MYYLIM(1);%*.9542;
plot(TimesToPlot-FRAMEDELAY,ValuesToPlotRescaled,'LineWidth',3,'Color','b');

% Cosmetics
xlim([0,42]);
ylim(MYYLIM);

title('Constitutive promoter dynamics');
MW_makeplotlookbetter(20);
%}