% CRP.cAMP project
% ===
% MW 2015/03/27

% This file plots Philippe's data; it is an extension of plotting Noreen's
% data, since philippe's data is only available as raw Schnitzcells files.
% All parameters are thus extracted immediately from the Schnitz files.

% Plot some of Noreen's extr/intr noise data, but with different aim.
% Data is available at:
% F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal
% Data is extracted from:
% \\biofysicasrv\Users2\Walker\ExtrinsicNoise\Data_Collection\Data
% See for more info:
% \\biofysicasrv\Users2\Walker\ExtrinsicNoise\Data_Collection\List_Experiments_Used_and_Unsused_ExtrNoiseProject.xlsx
% (Note that 1/3 of cellcycle means N datapoints used, which can be seen
% from name mu parameter.)

%% Define datasets

% Note that CRP.cAMP regulates when there's no induction; i.e. basal 
% expression should show ~different behavior as induced.
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Maltose basal 2012-05-16 pos1.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Maltose basal 2012-07-26 pos4.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Acetate basal 2012-07-19 pos4.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\RDM basal 2011-10-19 pos1.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\RDM basal 2011-10-19 pos3.mat');


% The datasets used for full induction figure
myDataFiles = ...
 {'F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\acetate-2011-07-29-pos1crop-Schnitz.mat', ...  
 'F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\lactose-2012-01-13-pos5crop-Schnitz.mat', ...  
 'F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\succinate-2011-09-02-pos3crop-Schnitz.mat', ...  
 'F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\defined-rich-2011-05-12-pos2crop-Schnitz.mat' ...  
     };

myGrouping = [ 1,2,3,4 ...
   ];

% These datasets have different fields for the parameters we're after
% (which indeed also means they might be determined slightly differently).
associatedFieldNames = ...
[ ...
 {'muP23_fitNew_cycCor','dY5_sum_dt_cycCor','Y6_mean_cycCor'}; ...
 {'muP19_fitNew_cycCor','dY5_sum_dt_cycCor','Y6_mean_cycCor'}; ...
 {'muP37_fitNew_cycCor','dY5_sum_dt_cycCor','Y6_mean_cycCor'}; ...
 {'muP15_fitNew_cycCor','dY5_sum_dt_cycCor','Y6_mean_cycCor'} ...
];

% I have extracted the following saved fitTimes from the Excel files in 
% the directories of the datasets. (See
% A_README_list-of-original-locations.txt).
load('F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\fitTimes.mat');
% Convert them to something processable (note this code works only when
% each sugar has only one dataset.
theFitTimes = cell2mat([fitTimes(find(strcmp(fitTimes(:,1),'acetate')),2);...
                fitTimes(find(strcmp(fitTimes(:,1),'lactose')),2);...
                fitTimes(find(strcmp(fitTimes(:,1),'succinate')),2);...
                fitTimes(find(strcmp(fitTimes(:,1),'DRM')),2)...
                ]);

legendDescriptions = {'ace','lac ','succ','DRM'};

some_colors; % loads some default color settings etc.

% TODO: check illum time!
 
%% Simple params
cyan    =  [0,1,1]/2;
yellow  =  [212/255, 170/255, 0];

%% First obtain all desired data

% Loop over desired datafiles
numberOfDataFiles = numel(myDataFiles);
myData = struct([]);
for i = 1:numberOfDataFiles
    % Load the datafile
    % clear mu_Third_cycCor, dY5_cycCor, fitydY, Y5_mean_cycCor;
    load(myDataFiles{i});
    
    % Prepare dummy p struct
    p={}, p.movieName = 'dummy';
    
    % Retrieve desired data
    % - mu (multiple points / individual)
    % - enzyme production rate (only Y; since has same filterset)
    % - enzyme concentration (only Y; since has same filterset)
    % - fitted line mu & enzyme expression
    allTimes = [schnitzcells(:).time];

    % Retrieve ALL growth rates
    retrievedData = DJK_get_schnitzData(p, schnitzcells, 'time','dataFields',{associatedFieldNames{i,1}},'fitTime',theFitTimes(i,:) );
    mu_values_all = [retrievedData.(associatedFieldNames{i,1})];    
   
    % Retrieve ALL production rates - not sure which field to use..!!
    retrievedData = DJK_get_schnitzData(p, schnitzcells, 'time','dataFields',{associatedFieldNames{i,2}},'fitTime',theFitTimes(i,:) );
    production_rates_all = [retrievedData.(associatedFieldNames{i,2})];    
    
    % Get indices of entries which correspond to fluor measurement, so we
    % can select above data for only those frames where there was an actual
    % fluor measurement.
    fluorIdx = find(~isnan(production_rates_all));
    % Update retrieved data using this information
    selected_growth_rates = mu_values_all(fluorIdx);
    selected_production_rates = production_rates_all(fluorIdx); 
    
    %NOW HAVE production rate: dY5_sum_dt_cycCor;
    %NOW HAVE growth rate: muP23_fitNew_cycCor      
        
    meanGrowthRate = mean(selected_growth_rates);
    
    % fit data and plot
    toFitGrowthRateRange = [prctile(selected_growth_rates,20),prctile(selected_growth_rates,80)];
    pdY = polyfit(selected_growth_rates',selected_production_rates',1);
    fitydY = pdY(1)*toFitGrowthRateRange + pdY(2);
    %NOW HAVE production fitline: fitydY
        
    % Retrieve data Y5_mean_cycCor
    retrievedData = DJK_get_schnitzData(p, schnitzcells, 'time','dataFields',{associatedFieldNames{i,3}},'fitTime',theFitTimes(i,:) );
    concentrations_all = [retrievedData.(associatedFieldNames{i,3})];
    % Select the ones w. fluorcolor
    selected_concentrations = concentrations_all(fluorIdx); 
    % NOW HAVE concentration: Y6_mean_cycCor        
    
    % fitline concentration
    pY = polyfit(selected_growth_rates',selected_concentrations',1);
    fityY = pY(1)*toFitGrowthRateRange + pY(2);
    % NOW HAVE concentration fitline: fityY
    
    % Save data to convenient plotting struct
    myData(i).descrip = myDataFiles(i);
    myData(i).selected_growth_rates = selected_growth_rates;
    myData(i).selected_production_rates = selected_production_rates;
    myData(i).selected_concentrations = selected_concentrations;
    myData(i).fullGrowthRateRange = toFitGrowthRateRange;
    myData(i).fitydY = fitydY;
    myData(i).fityY = fityY;
    myData(i).averagemu = meanGrowthRate;    
    myData(i).averagedY = mean(selected_production_rates);
    myData(i).averageY = mean(selected_concentrations);
    
    % Tell user what we're doing
    disp(['==== Loading data ' num2str(i) ' of ' num2str(numberOfDataFiles) ' done. ====']);
end


%% Then plot all this data

% plotting production rate plot
hFig = figure(1), clf; hold on;
offset=100; width1=1000; height1=500;
set(hFig, 'Position', [offset offset width1 height1])
subplot(1,2,1), hold on;

for i = 1:numberOfDataFiles    
    % cloud
    plot(myData(i).selected_growth_rates,myData(i).selected_production_rates,'o','Color',preferredcolors(myGrouping(i)+1,:));
end

for i = 1:numberOfDataFiles    
    % fit
    plot(myData(i).fullGrowthRateRange,myData(i).fitydY,'-','Color','k','LineWidth',3)
end

% average point (used for legend too)
legendLines = []; previous = 0;
for i = 1:numberOfDataFiles        
    lineH = plot(myData(i).averagemu,myData(i).averagedY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15)
if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','best');


xlabel('Growth rate (dbl/hr)');
ylabel('Producation rate');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

allMu=[myData(:).selected_growth_rates];
xlim([0, max(allMu)]);

% plotting concentration plot
% =================================
hFig = figure(1); hold on;
%offset=100; width2=500; height2=500;
%set(hFig, 'Position', [offset*2+width1 offset width2 height2])
subplot(1,2,2), hold on;

% cloud
for i = 1:numberOfDataFiles        
    plot(myData(i).selected_growth_rates,myData(i).selected_concentrations,'o','Color',preferredcolors(myGrouping(i)+1,:));    
end

% fit
for i = 1:numberOfDataFiles        
    plot(myData(i).fullGrowthRateRange,myData(i).fityY,'-','Color','k','LineWidth',3)
end

% average point
%legendLines = []; previous = 0;
for i = 1:numberOfDataFiles        
    lineH = plot(myData(i).averagemu,myData(i).averageY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15)
    %if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
%legend( legendLines, {'ace','lac ','mal','RDM'},'Location','best');

xlabel('Growth rate (dbl/hr)');
ylabel('Concentration');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

allMu=[myData(:).selected_growth_rates];
xlim([0, max(allMu)]);







