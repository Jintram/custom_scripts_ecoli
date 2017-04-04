
NOSAVEPLEASE=1;
YLIMbirthlife=[0,200];

myPlotting.style = 'CBmanuscript';

%% dataset I
WHATDATA = 'tetracycline';
LENGTHFIELD = 'length_fitNew';

RUNSECTIONSFILADIV = 'loadData';
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'lifeTimeVsBirthLength2';
script20160429_filamentRecoveryDivisionRatioss

combinedDynamicsData.(WHATDATA).binCentersDynData = binCentersDynData;
combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData = meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData = stdValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).CVForBinsDynData = stdValuesForBinsDynData./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).VarByMeanForBinsDynData = stdValuesForBinsDynData.^2./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).selectedXdata = selectedXdata; % selected for after switchtime
combinedDynamicsData.(WHATDATA).selectedYdata = selectedYdata;

%% dataset II
WHATDATA = 'temperature';
LENGTHFIELD = 'length_skeleton';

RUNSECTIONSFILADIV = 'loadData';
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'lifeTimeVsBirthLength2';
script20160429_filamentRecoveryDivisionRatioss

combinedDynamicsData.(WHATDATA).binCentersDynData = binCentersDynData;
combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData = meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData = stdValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).CVForBinsDynData = stdValuesForBinsDynData./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).VarByMeanForBinsDynData = stdValuesForBinsDynData.^2./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).selectedXdata = selectedXdata; % selected for after switchtime
combinedDynamicsData.(WHATDATA).selectedYdata = selectedYdata;

%% dataset III
WHATDATA = 'sulA';
LENGTHFIELD = 'length_skeleton';

RUNSECTIONSFILADIV = 'loadData';
script20160429_filamentRecoveryDivisionRatioss
RUNSECTIONSFILADIV = 'lifeTimeVsBirthLength2';
script20160429_filamentRecoveryDivisionRatioss

combinedDynamicsData.(WHATDATA).binCentersDynData = binCentersDynData;
combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData = meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData = stdValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).CVForBinsDynData = stdValuesForBinsDynData./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).VarByMeanForBinsDynData = stdValuesForBinsDynData.^2./meanValuesForBinsDynData;
combinedDynamicsData.(WHATDATA).selectedXdata = selectedXdata; % selected for after switchtime
combinedDynamicsData.(WHATDATA).selectedYdata = selectedYdata;

%% Plotting together

h=figure; clf; hold on; 

myColors = linspecer(3);

DATATOPLOT = {'tetracycline', 'temperature', 'sulA'}
%DATATOPLOT = {'tetracycline'}

legendLines=[];
% Plot large circle + error bar if multiple datapoints
for dataSetindex = 1:numel(DATATOPLOT)
    
    % Parameter for dataset selection
    WHATDATA=DATATOPLOT{dataSetindex}
    % And color
    theColor = myColors(dataSetindex,:);
    
    % Select data that have multiple points
    pointsWithMultipleDataIndices = combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData>0
    
    % Plots those
    xdata   = combinedDynamicsData.(WHATDATA).binCentersDynData(pointsWithMultipleDataIndices);
    ydata   = combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData(pointsWithMultipleDataIndices);
    errdata = combinedDynamicsData.(WHATDATA).stdValuesForBinsDynData(pointsWithMultipleDataIndices);
    errorbar(xdata,ydata,errdata,...
            '-o','Color',theColor,'MarkerSize',5,'LineWidth',2);%,'LineWidth',1)
    %plot(combinedDynamicsData.(WHATDATA).binCentersDynData,combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData,'.','Color',theColor,'LineWidth',2);
    combinedDynamicsData.(WHATDATA).binCentersDynData(1)
end
% single data on top
for dataSetindex = 1:numel(DATATOPLOT)
    
    WHATDATA=DATATOPLOT{dataSetindex};
    theColor = myColors(dataSetindex,:);    
    
    %legendLines(end+1)=plot(combinedDynamicsData.(WHATDATA).binCentersDynData,combinedDynamicsData.(WHATDATA).meanValuesForBinsDynData,'x','Color',theColor,'LineWidth',3,'MarkerSize',10);    
    
    selectedXdata = combinedDynamicsData.(WHATDATA).selectedXdata;
    selectedYdata = combinedDynamicsData.(WHATDATA).selectedYdata;
    %scatter(selectedXdata,selectedYdata,'.','Color',theColor); % selected for after switchtime

end

ylim([0,100]);
xlabel('Birth size (\mum)');
ylabel('Interdivision time (min)');
MW_makeplotlookbetter(20);

legend(legendLines,{'Tetracycline', 'Temperature', 'SulA'});

hFigSIIC = h;
