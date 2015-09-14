

%% Parameter settings
USERSETTINGS.myRootDir='U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\';
USERSETTINGS.myDateDir='2015_09_01_exp1\';
USERSETTINGS.datafile='2015_09_01_dataCRPexperiment.xls';
USERSETTINGS.customSuffix = '_GFP';
USERSETTINGS.ODorFluor = 2;


myConcentrationValues = [  ...
        0.001 ,...
        0.0003 ,...
        0.00009 ,...
        0.000027 ,...
        0.0000081 ,...
        0
        ];

myXlim = [1e-6,max(myConcentrationValues)*10];

% Which points to take into account for fit, in relation to
% myConcentrationValues defined above
STARTAT = 1;
ENDIGNORE = 0;


%% Call script for first batch
DONTSAVE=1;
USERSETTINGS.wellNamesToPlot = ...
    {'str38A', 'str38B', 'str38C', 'str38D', 'str38E', 'str38F'};
ExtractFitPlateReaderData_General_Part2_Fluor

%%
figure(88), clf, hold on
legendLinesHandles = [];

% y=zero line
plot(myXlim, [0, 0], 'k-','LineWidth',2)

for i = 1:numel(output)
    l=plot(   ones(1,numel(output(i).manualMuValues)) * myConcentrationValues(i), ...
                output(i).manualMuValues,'o')
    set(l, 'MarkerSize', 10,'LineWidth',3,'Color','b');   
end
legendLinesHandles(end+1) = l;
% value at zero
plot(myXlim,[1,1]*mean(output(end).manualMuValues),'--b','LineWidth',3);

% Some settings for the plot
set(gca,'xscale','log');

MW_makeplotlookbetter(20);

ylim([-.1,1]);
xlim(myXlim);

xlabel('cAMP concentration [M]');
ylabel('fitted growth rate [dbl/hr]');

%% now retrieve the cAMP extracell. data

DONTSAVE=1;
USERSETTINGS.wellNamesToPlot = ...
    {'str39A', 'str39B', 'str39C', 'str39D', 'str39E', 'str39F'};
ExtractFitPlateReaderData_General_Part2_Fluor

%% and plot it
figure(88), hold on;
for i = 1:numel(output)
    l=plot(   ones(1,numel(output(i).manualMuValues)) * myConcentrationValues(i), ...
                output(i).manualMuValues,'s')
    set(l, 'MarkerSize', 10,'LineWidth',3,'Color', 'r');
end
legendLinesHandles(end+1) = l;
% value at zero
plot(myXlim,[1,1]*mean(output(end).manualMuValues),'--r','LineWidth',3);

% Make a fit
%{
meanValues = [];
disp('Fitting, ignoring last point.')
for i = STARTAT:(numel(output)-ENDIGNORE) % last point, [cAMP] zero
    meanValues(end+1) = mean(output(i).manualMuValues);    
end
toFitX = myConcentrationValues(STARTAT:end-ENDIGNORE)/log(10);
p=polyfit(toFitX,meanValues,2)

fittedXs=linspace(myXlim(1),myXlim(2),100);
logFittedXs=log(fittedXs)/log(10);
fittedYs=p(1)*logFittedXs.^2 + p(2).*logFittedXs + p(3);
%figure
plot(fittedXs,fittedYs,'-r','LineWidth',2);
set(gca,'xscale','log');

% Gaussian fit
f = fit(toFitX',meanValues','gauss2')
%}
%% 

legend(legendLinesHandles, {'838..','839..'});




