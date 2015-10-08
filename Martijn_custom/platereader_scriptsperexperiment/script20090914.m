


%% Parameter settings
USERSETTINGS.myRootDir='U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\';
USERSETTINGS.myDateDir='2015_09_01_exp1\';
USERSETTINGS.datafile='2015_09_01_dataCRPexperiment.xls';
USERSETTINGS.customSuffix = '_GFP';
USERSETTINGS.ODorFluor = 2;


USERSETTINGS.myConcentrationValues = ... % 0.0031 M start, % 3.1622X dilution series
        [(94/2994*.1) ./ (3.1622 .^ [0 1 2 3 4 5 6]) 0]

USERSETTINGS.myXlim = [1e-6,max(USERSETTINGS.myConcentrationValues)*10];

% Which points to take into account for fit, in relation to
% myConcentrationValues defined above
USERSETTINGS.STARTAT = 1;
USERSETTINGS.ENDIGNORE = 0;


%% Call script for first batch
myhandle=figure(88), clf, hold on;
legendLinesHandles = [];

% s70 data WT strain
DONTSAVE=1;
USERSETTINGS.wellNamesToPlot = ...
{'str841A','str841B','str841C','str841D','str841E','str841F','str841G','str841H'};
ExtractFitPlateReaderData_General_Part2_Fluor
[allLinehandles, legendLinesHandles] = platereaderextraplot(USERSETTINGS,output,myhandle,legendLinesHandles,[.6 .6 1])
for l = allLinehandles
    set(l, 'Marker', 's', 'MarkerSize', 10,'LineWidth',2);
end

% s70 data delta strain
DONTSAVE=1;
USERSETTINGS.wellNamesToPlot = ...
{'str894A','str894B','str894C','str894D','str894E','str894F','str894G','str894H'};
ExtractFitPlateReaderData_General_Part2_Fluor
[allLinehandles, legendLinesHandles] = platereaderextraplot(USERSETTINGS,output,myhandle,legendLinesHandles,[1 .6 .6])
for l = allLinehandles
    set(l, 'Marker', 's', 'MarkerSize', 10,'LineWidth',2);
end

% CRP data WT strain
DONTSAVE=1;
USERSETTINGS.wellNamesToPlot = ...
{'str842A','str842B','str842C','str842D','str842E','str842F','str842G','str842H'};
ExtractFitPlateReaderData_General_Part2_Fluor
[allLinehandles, legendLinesHandles] = platereaderextraplot(USERSETTINGS,output,myhandle,legendLinesHandles,'b')
% CRP data delta strain 
DONTSAVE=1;
USERSETTINGS.wellNamesToPlot = ...
{'str893A','str893B','str893C','str893D','str893E','str893F','str893G','str893H'};
ExtractFitPlateReaderData_General_Part2_Fluor
[allLinehandles, legendLinesHandles] = platereaderextraplot(USERSETTINGS,output,myhandle,legendLinesHandles,'r')

%***

legend(legendLinesHandles,{'s70 WT','s70 delta','CRP WT','CRP delta'},'Location','northeastoutside')

%% And plot






