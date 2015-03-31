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

myDataFiles = ...
 {'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Acetate full 2012-06-15 pos4.mat', ...
  ... % 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Acetate full 2012-06-20 pos2.mat', ... % Sanity check failed!
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-01-27 pos1.mat', ...
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-03-02 pos6.mat', ...
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Lactose full 2012-04-02 pos4.mat', ...  
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full 2012-04-24 pos5.mat', ...
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full 2012-05-08 pos2.mat', ...
  ... % 'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\Maltose full highIllum incl more bimodal 2012-05-24 pos5.mat', ...% Sanity check failed!
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\RDM full 2011-10-28 pos1.mat', ...
  'F:\X_Other_datasets\CRPcAMP\NWdata-lacfullinduced\RDM full 2011-10-28 pos2.mat' ...  
     };

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
    myData(i).mu_Third_cycCor = mu_Third_cycCor;
    myData(i).dY5_cycCor = dY5_cycCor;
    myData(i).Y5_mean_cycCor = Y5_mean_cycCor;
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
hFig = figure(1), clf;
offset=100; width=1000; height=500;
set(hFig, 'Position', [offset offset width height])
subplot(1,2,1), hold on;

for i = 1:numberOfDataFiles    
    % cloud
    plot(myData(i).mu_Third_cycCor,myData(i).dY5_cycCor,'o','Color',preferredcolors(myGrouping(i)+1,:));
end

for i = 1:numberOfDataFiles    
    % fit
    plot(myData(i).fullGrowthRateRange,myData(i).fitydY,'-','Color','k','LineWidth',3)
end

for i = 1:numberOfDataFiles    
    % average point
    plot(myData(i).averagemu,myData(i).averagedY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15)
end    

xlabel('Growth rate (dbl/hr)');
ylabel('Producation rate');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

allMu=[myData(:).mu_Third_cycCor];
xlim([0, max(allMu)]);

%% plotting concentration plot
hFig = figure(1);
offset=100; width=1000; height=500;
set(hFig, 'Position', [offset offset width height])
subplot(1,2,2), hold on;

% cloud
for i = 1:numberOfDataFiles        
    plot(myData(i).mu_Third_cycCor,myData(i).Y5_mean_cycCor,'o','Color',preferredcolors(myGrouping(i)+1,:));    
end

% fit
for i = 1:numberOfDataFiles        
    plot(myData(i).fullGrowthRateRange,myData(i).fityY,'-','Color','k','LineWidth',3)
end

% average point
legendLines = [];
for i = 1:numberOfDataFiles        
    plot(myData(i).averagemu,myData(i).averageY,'o','MarkerFaceColor',preferredcolors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',15)
    if myGrouping(i) ~= previous, legendLines = [legendLines i]; end; previous = myGrouping(i);
end    
%legend( legendLines+numberOfDataFiles*2,{'ace','lac ','mal','RDM'},'Location','best');%,'Location','southeast')
legend( legendLines+numberOfDataFiles*2,{'ace','lac ','mal','RDM'},'Location','southeast');%,'Location','southeast')

xlabel('Growth rate (dbl/hr)');
ylabel('Concentration');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

allMu=[myData(:).mu_Third_cycCor];
xlim([0, max(allMu)]);







