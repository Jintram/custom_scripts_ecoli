% CRP.cAMP project
% ===
% MW 2015/03/27

% This file plots Philippe's data; it is an extension of plotting Noreen's
% data, since philippe's data is only available as raw Schnitzcells files.
% All parameters are thus extracted immediately from the Schnitz files.

% I copied the datasets that seemed to correspond to Philippe's measurement
% to F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\
% Where also a list of all experiments can be found.
%
% Kiviet 2014 et al. doesn't mention that these assays were performed under
% IPTG (in fact, it states IPTG is only used 'when indicated'). 

%% Define datasets

% The datasets used for full induction figure based Kiviet2014
myDataFiles = ...
 {'F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\acetate-2011-07-29-pos1crop-Schnitz.mat', ...  
 'F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\lactose-2012-01-13-pos5crop-Schnitz.mat', ...  
 ... 'F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\succinate-2011-09-02-pos2crop-Schnitz.mat', ...  
 'F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\succinate-2011-09-02-pos3crop-Schnitz.mat', ...  
 'F:\X_Other_datasets\CRPcAMP\PNdata-lacfullinduced\defined-rich-2011-05-12-pos2crop-Schnitz.mat' ...  
     }

myGrouping = [ 1,2,3,4 ...
   ];

% These datasets have different fields for the parameters we're after
% (which indeed also means they might be determined slightly differently).
associatedFieldNames = ...
[ ...
%{
{'muP23_fitNew_cycCor','dY5_sum_dt_cycCor','Y5_mean'}; ...
 {'muP19_fitNew_cycCor','dY5_sum_dt_cycCor','Y5_mean'}; ...
 {'muP37_fitNew_cycCor','dY5_sum_dt_cycCor','Y5_mean'}; ...
 {'muP15_fitNew_cycCor','dY5_sum_dt_cycCor','Y5_mean'} ...
%}
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
                fitTimes(find(strcmp(fitTimes(:,1),'succinatepos3')),2);...
                fitTimes(find(strcmp(fitTimes(:,1),'DRM')),2)...
                ]);

legendDescriptions = {'Acetate','Lactose ','Succinate','Rich defined'};

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

    % Retrieving data I
    % ===
    
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
    
    % Fitting growth rate data 
    % ===
    
    % E(mu), 2nd order  
    toFitGrowthRateRange = [prctile(selected_growth_rates,5),prctile(selected_growth_rates,95)];    
    N = 100;
    fitlineratesXmus = [toFitGrowthRateRange(1):(toFitGrowthRateRange(2)-toFitGrowthRateRange(1))/N:toFitGrowthRateRange(2)];
    coefficientsFitLinesMuAsX = polyfit(selected_growth_rates',selected_production_rates',2);
    fitlineRatesYE = polyval(coefficientsFitLinesMuAsX,fitlineratesXmus); % coefficientsFitLinesMuAsX(1)*fitlineratesXmus.^2 + coefficientsFitLinesMuAsX(2)*fitlineratesXmus + coefficientsFitLinesMuAsX(3);
    %NOW HAVE linear production fitline: fitydY
    
    % mu(E), 2nd order polyfit   
    toFitEnzymeRateRange = [prctile(selected_production_rates,5),prctile(selected_production_rates,95)];
    N = 100;
    dEnzyme = (toFitEnzymeRateRange(2)-toFitEnzymeRateRange(1))/N;
    fitlineratesXenzyme = [toFitEnzymeRateRange(1):dEnzyme:toFitEnzymeRateRange(2)];
    coefficientsFitLinesEnzymeAsX = polyfit(selected_production_rates',selected_growth_rates',2);
    fitlineRatesYmu = polyval(coefficientsFitLinesEnzymeAsX,fitlineratesXenzyme); % (1)*fitlineratesXenzyme.^2 + coefficientsFitLinesEnzymeAsX(2)*fitlineratesXenzyme + coefficientsFitLinesEnzymeAsX(3);
    % NOW HAVE polynomial fit mu(E) = a E^2 + b E + c
    
    % Retrieving data II
    % ===
        
    % Retrieve data Y5_mean_cycCor
    retrievedData = DJK_get_schnitzData(p, schnitzcells, 'time','dataFields',{associatedFieldNames{i,3}},'fitTime',theFitTimes(i,:) );
    concentrations_all = [retrievedData.(associatedFieldNames{i,3})];
    % Select the ones w. fluorcolor
    selected_concentrations = concentrations_all(fluorIdx); 
    % NOW HAVE concentration: Y6_mean_cycCor        
    
    % fitting concentrations
    % ===
    
    % 2nd order polyfit, E as function mu: E(mu)
    N = 100;
    dMu = (toFitGrowthRateRange(2)-toFitGrowthRateRange(1))/N;
    fitlinesConcXmu = [toFitGrowthRateRange(1):dMu:toFitGrowthRateRange(2)];
    coefficientsFitLinesConcMuAsX = polyfit(selected_growth_rates',selected_concentrations',2);% TODO make this arbitrary order
    fitlineConcYE = polyval(coefficientsFitLinesConcMuAsX,fitlinesConcXmu); %coefficientsFitLinesConcMuAsX(1)*fitlinesConcXmu.^2 + coefficientsFitLinesConcMuAsX(2)*fitlinesConcXmu + coefficientsFitLinesConcMuAsX(3);
    % NOW HAVE concentration fitline: fityY

    % 2nd order polyfit, mu as function E: mu(E)
    toFitEnzymeConcentrationRange = [prctile(selected_concentrations,5),prctile(selected_concentrations,95)];
    N = 100;
    dE = (toFitEnzymeConcentrationRange(2)-toFitEnzymeConcentrationRange(1))/N;
    fitlinesConcXenzymes = [toFitEnzymeConcentrationRange(1):dE:toFitEnzymeConcentrationRange(2)];
    coefficientsFitLinesConcEnzymeAsX = polyfit(selected_concentrations',selected_growth_rates',2);
    fitlineConcYmu = polyval(coefficientsFitLinesConcEnzymeAsX,fitlinesConcXenzymes); % (1)*fitlinesConcXenzymes.^2 + coefficientsFitLinesConcEnzymeAsX(2)*fitlinesConcXenzymes + coefficientsFitLinesConcEnzymeAsX(3);
    % NOW HAVE polynomial fit mu(E) = a E^2 + b E + c    
    
    % Save data to convenient plotting struct
    % ===
    myData(i).descrip = myDataFiles(i);
    myData(i).selected_growth_rates = selected_growth_rates;
    myData(i).selected_production_rates = selected_production_rates;
    myData(i).selected_concentrations = selected_concentrations;

    myData(i).averagemu = meanGrowthRate;    
    myData(i).averagedY = mean(selected_production_rates);
    myData(i).averageY = mean(selected_concentrations);
        
    % Rate fits
    % ===  
    % E(mu)
    myData(i).fitlineratesXmus = fitlineratesXmus;
    myData(i).coefficientsFitLinesMuAsX = coefficientsFitLinesMuAsX;
    myData(i).fitlineRatesYE = fitlineRatesYE;
    % mu(E)
    myData(i).fitlineratesXenzyme = fitlineratesXenzyme;
    myData(i).coefficientsFitLinesEnzymeAsX = coefficientsFitLinesEnzymeAsX;
    myData(i).fitlineRatesYmu = fitlineRatesYmu;    
    
    % Concentration fits
    % ===
    % E(mu)
    myData(i).fitlinesConcXmu = fitlinesConcXmu;
    myData(i).coefficientsFitLinesConcMuAsX = coefficientsFitLinesConcMuAsX;
    myData(i).fitlineConcYE = fitlineConcYE;
    % mu(E)
    myData(i).fitlinesConcXenzymes = fitlinesConcXenzymes;
    myData(i).coefficientsFitLinesConcEnzymeAsX = coefficientsFitLinesConcEnzymeAsX;
    myData(i).fitlineConcYmu = fitlineConcYmu;
    
    % Tell user what we're doing
    disp(['==== Loading data ' num2str(i) ' of ' num2str(numberOfDataFiles) ' done. ====']);
end


%% Then plot all this data

% plotting production rate plot
hFig = figure(1), clf; hold on;
offset=100; width1=1000; height1=500;
set(hFig, 'Position', [offset offset width1 height1]);
subplot(1,2,1), hold on;

% cloud
for i = 1:numberOfDataFiles        
    plot(myData(i).selected_growth_rates,myData(i).selected_production_rates,'o','Color',preferredcolors(myGrouping(i)+1,:));
end

% fit
for i = 1:numberOfDataFiles
    % E(mu)
    plot(myData(i).fitlineratesXmus,myData(i).fitlineRatesYE,'-','Color','k','LineWidth',3)
    % mu(E)
    plot(myData(i).fitlineRatesYmu,myData(i).fitlineratesXenzyme,'-','Color','r','LineWidth',3)        
end

% average point (used for legend too)
legendLines = []; previous = 0;
for i = 1:numberOfDataFiles        
    lineH = plot(myData(i).averagemu,myData(i).averagedY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15);
if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','northeast');


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
    plot(myData(i).selected_growth_rates,myData(i).selected_concentrations,'.','Color',preferredcolors(myGrouping(i)+1,:),'MarkerSize',3);    
end

% fit
for i = 1:numberOfDataFiles        
    % E(mu)
    plot(myData(i).fitlinesConcXmu,myData(i).fitlineConcYE,'-','Color','k','LineWidth',3)
    % mu(E)
    plot(myData(i).fitlineConcYmu,myData(i).fitlinesConcXenzymes,'-','Color','r','LineWidth',3)    
end


% average point
legendLines = []; previous = 0;
for i = 1:numberOfDataFiles        
    lineH = plot(myData(i).averagemu,myData(i).averageY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15);
    if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','northeast');

xlabel('Growth rate (dbl/hr)');
ylabel('Concentration');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

allMu=[myData(:).selected_growth_rates];
xlim([0, max(allMu)]);


%% Kernel plot
% ==========

hFig = figure(2); clf;
offset=100; width1=1000; height1=500;
set(hFig, 'Position', [offset offset width1 height1]);
subplot(1,2,1); hold on;

% Growth rate data
% ===
for i = 1:numberOfDataFiles
    data = [myData(i).selected_growth_rates', myData(i).selected_production_rates'];
    [bandwidth,density,X,Y] = kde2d(data);    
    [C, l1] = contour3(X,Y,density,2,'-k');
    set(l1, 'LineWidth', 1);
    plot(data(:,1),data(:,2),'.','Color',preferredcolors(myGrouping(i)+1,:),'MarkerSize',3);
end

% average point (used for legend too)
legendLines = []; previous = 0;
for i = 1:numberOfDataFiles        
    lineH = plot(myData(i).averagemu,myData(i).averagedY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);
if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','northeast');

xlabel('Growth rate (dbl/hr)');
ylabel('Production rate (a.u./min)');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);
ylim([-750, 2000])

%set(gca,'YScale','log');

% Concentration data
% ===
subplot(1,2,2); hold on;
for i = 1:numberOfDataFiles
    data = [myData(i).selected_growth_rates', myData(i).selected_concentrations'];
    [bandwidth,density,X,Y] = kde2d(data);    
    [C, l1] = contour3(X,Y,density,2,'-k'); 
    set(l1, 'LineWidth', 2);
    plot(data(:,1),data(:,2),'.','Color',preferredcolors(myGrouping(i)+1,:),'MarkerSize',3);
end

% averages
legendLines = []; previous = 0;
for i = 1:numberOfDataFiles        
    lineH = plot(myData(i).averagemu,myData(i).averageY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);
    if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','northeast');

%set(gca,'YScale','log');

ylim([0, 400])
xlabel('Growth rate (dbl/hr)');
ylabel('Concentration (a.u.)');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);






