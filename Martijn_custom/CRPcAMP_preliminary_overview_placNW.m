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

%% Define datasets

% Note that CRP.cAMP regulates when there's no induction; i.e. basal 
% expression should show ~different behavior as induced.
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Maltose basal 2012-05-16 pos1.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Maltose basal 2012-07-26 pos4.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Acetate basal 2012-07-19 pos4.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\RDM basal 2011-10-19 pos1.mat');
%load('F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\RDM basal 2011-10-19 pos3.mat');

% Illumination times!
% YFP, CFP



% The datasets used for full induction figure
myDataFiles = ...
 {'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Acetate full 2012-06-15 pos4.mat', ... % illum 20,40
  ... % 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Acetate full 2012-06-20 pos2.mat', ... % Sanity check failed!
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-01-27 pos1.mat', ... % illum 25, 100
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-03-02 pos6.mat', ... % illum 25, 100
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-04-02 pos4.mat', ... % illum 25, 100  
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full 2012-04-24 pos5.mat', ... % illum 25, 100
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full 2012-05-08 pos2.mat', ... % illum 25 ,30
  ... % 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full highIllum incl more bimodal 2012-05-24 pos5.mat', ...% Sanity check failed!
   'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\RDM full 2011-10-28 pos1.mat', ... % illum 25, 100
   'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\RDM full 2011-10-28 pos2.mat' ...  % illum 25, 100
     };
 % I used YFP, which always has ~25ms illumination time..

myGrouping = [ ...
      1, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Acetate full 2012-06-15 pos4.mat', ...
      ... % 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Acetate full 2012-06-20 pos2.mat', ... % Sanity check failed!
      2, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-01-27 pos1.mat', ...
      2, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-03-02 pos6.mat', ...
      3, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-04-02 pos4.mat', ...  
      3, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full 2012-04-24 pos5.mat', ...
      3, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full 2012-05-08 pos2.mat', ...
      ... % 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full highIllum incl more bimodal 2012-05-24 pos5.mat', ...% Sanity check failed!
      4, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\RDM full 2011-10-28 pos1.mat', ...
      4 ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\RDM full 2011-10-28 pos2.mat' ... 
    ];

legendDescriptions = {'Acetate','Lactose','Maltose','Rich defined'};


%{
% Dataset basal expression (only induced should respond to CRP.cAMP
% directly)
myDataFiles = ...
 {...'F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Acetate basal 2012-07-19 pos4.mat', ...
 'F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Maltose basal 2012-05-16 pos1.mat', ...
 'F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\Maltose basal 2012-07-26 pos4.mat', ...
 'F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\RDM basal 2011-10-19 pos1.mat', ...
 'F:\X_Other_datasets\CRPcAMP\NWdata-lacbasal\RDM basal 2011-10-19 pos3.mat', ...
     };

myGrouping = [ ...
      ...1, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Acetate full 2012-06-15 pos4.mat', ...
      2, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-01-27 pos1.mat', ...
      2, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-03-02 pos6.mat', ...
      3, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-04-02 pos4.mat', ...  
      3, ... 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full 2012-04-24 pos5.mat', ...      
    ];

legendDescriptions = {'mal','RDM'};
%}

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
    
    % Retrieve desired data
    % - mu (multiple points / individual)
    % - enzyme production rate (only Y; since has same filterset)
    % - enzyme concentration (only Y; since has same filterset)
    % - fitted line mu & enzyme expression
    %NOW HAVE growth rate: mu_Third_cycCor;
    %NOW HAVE production rate: dY5_cycCor;
    meanGrowthRate = mean(mu_Third_cycCor);
    
    % fit data and plot
    fullGrowthRateRange = [0, meanGrowthRate];
    pdY = polyfit(mu_Third_cycCor',dY5_cycCor',1);
    fitydY = pdY(1)*fullGrowthRateRange + pdY(2);
    %NOW HAVE production fitline: fitydY
    
    % Extract the concentration data using Daan's function
    % ===
    
    % Prepare dummy p struct
    p={}, p.movieName = 'dummy';
    % Retrieve data Y5_mean_cycCor
    retrievedData = DJK_get_schnitzData(p, schnitzcells_rm, 'time_atY','dataFields',{'Y5_mean_cycCor'},'fitTime',fitTime);
    Y5_mean_cycCor = [retrievedData.Y5_mean_cycCor];
    % NOW HAVE concentration: Y5_mean_cycCor
    
    % Sanity check, check if I retrieve same data as Noreen.
    retrievedData = DJK_get_schnitzData(p, schnitzcells_rm, 'time_atY','dataFields',{'dC5_cycCor'},'fitTime',fitTime);
    mydC5_cycCor = [retrievedData.dC5_cycCor];
    if numel(mydC5_cycCor) ~= numel(dC5_cycCor)
        error(['Sanity check failed in #' num2str(i) '! Unequal # entries.']);
    else
        sanityCheck = mydC5_cycCor - dC5_cycCor;    
        if ~isempty(sanityCheck(find(sanityCheck>0)))
            error(['Sanity check failed in #' num2str(i) '! sanityCheck chould be empty']);
        end
        % check this also
    end
    
    % fitline concentration
    pY = polyfit(mu_Third_cycCor',Y5_mean_cycCor',1);
    fityY = pY(1)*fullGrowthRateRange + pY(2);
    % NOW HAVE concentration fitline: fityY
    
    % Save data to convenient plotting struct
    myData(i).descrip = descrip;
    myData(i).selected_growth_rates = mu_Third_cycCor;
    myData(i).selected_production_rates = dY5_cycCor;
    myData(i).selected_concentrations = Y5_mean_cycCor;
    
    myData(i).fullGrowthRateRange = fullGrowthRateRange;    
    myData(i).fitydY = fitydY;
    myData(i).fityY = fityY;
    myData(i).averagemu = meanGrowthRate;    
    myData(i).averagedY = mean(dY5_cycCor);
    myData(i).averageY = mean(Y5_mean_cycCor);
    
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


%% Kernel plot
% ==========

hFig = figure(2); clf;
offset=100; width1=1000; height1=500;
set(hFig, 'Position', [offset offset width1 height1]);
subplot(1,2,1); hold on;

% Growth rate data
% ===
for i = 1:(numberOfDataFiles-2)
    data = [myData(i).selected_growth_rates', myData(i).selected_production_rates'];
    [bandwidth,density,X,Y] = kde2d(data);    
    [C, l1] = contour3(X,Y,density,1,'-k');
    set(l1, 'LineWidth', 2,'Color', preferredcolors(myGrouping(i)+1,:)/2);
    plot(data(:,1),data(:,2),'.','Color',preferredcolors(myGrouping(i)+1,:),'MarkerSize',3);
end

% average point (used for legend too)
legendLines = []; previous = 0;
for i = 1:(numberOfDataFiles-2)
    lineH = plot(myData(i).averagemu,myData(i).averagedY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);
if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','northeast');

xlabel('Growth rate (dbl/hr)');
ylabel('Production rate (a.u./min)');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);
ylim([-750, 10000])
xlim([0, 1]);

%set(gca,'YScale','log');

% Concentration data
% ===
subplot(1,2,2); hold on;
for i = 1:(numberOfDataFiles-2)
    data = [myData(i).selected_growth_rates', myData(i).selected_concentrations'];
    [bandwidth,density,X,Y] = kde2d(data);    
    [C, l1] = contour3(X,Y,density,1,'-'); 
    set(l1, 'LineWidth', 2,'Color', preferredcolors(myGrouping(i)+1,:)/2);
    plot(data(:,1),data(:,2),'.','Color',preferredcolors(myGrouping(i)+1,:),'MarkerSize',3);
end

% averages
legendLines = []; previous = 0;
for i = 1:(numberOfDataFiles-2)        
    lineH = plot(myData(i).averagemu,myData(i).averageY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);
    if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','northeast');

%set(gca,'YScale','log');

%ylim([0, 400])
xlim([0, 1]);

xlabel('Growth rate (dbl/hr)');
ylabel('Concentration (a.u.)');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);






