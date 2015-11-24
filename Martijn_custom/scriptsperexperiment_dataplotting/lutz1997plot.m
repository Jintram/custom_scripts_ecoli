
%% Information from expressys website:
% http://www.expressys.com/main_applications.html
%- Expression levels of promoter PLtetO-1 show sigmoidal character between 1 and 20ng/ml aTc.
%- An expression plateau is reached above 20ng aTc.

%% Settings
MYMAXX = inf;

some_colors % script that loads colors

%% Lutz1997 plot extracted using plot digitizer (aka datathief?)

noinduction = [0,0.97355276];

% data
data=[...
0.9834712,1.1594992;...
1.9434516,1.720024;...
2.9646256,2.4280143;...
4.050459,5.4428782;...
4.931296,8.596791;...
6.005728,14.334181;...
6.752361,27.406288;...
8.083515,38.852993;...
9.693564,72.21745;...
11.426331,114.12942;...
12.444713,236.81851;...
16.568851,1486.9135;...
16.873125,2057.4722;...
20.845118,2545.8853;...
23.689047,2045.1888;...
24.90363,2607.6;...
36.08717,2324.5312;...
65.64396,2494.862;...
]

%% Select data & extract params!
mySelection = find(data(:,1)<MYMAXX);

allXData = data(:, 1)';
allYData = data(:, 2)';

selectedxData = data( mySelection , 1 )';
selectedyData = data( mySelection , 2 )';

highestYValue = max(allYData);

myYlimmin = min([allYData])/2; 
myYlimmax = max(allYData)*2;

%% plot


%%
hill = @(myParameters,xdata) ...
    myParameters(1) * xdata.^myParameters(2) ./ (xdata.^myParameters(2) + myParameters(3)) + myParameters(4);
    % myParameters(1) = ymax
    % myParameters(2) = n
    % myParameters(3) = K
    % myParameters(4) = offset

myParameters0=[highestYValue,1,1000,.1];
fittedParameters=nlinfit(selectedxData,selectedyData,hill,myParameters0)

%% Now to log transform (this is a bit redundant)
loghill = @(myParameters,xdata) ...
    log(myParameters(1) * xdata.^myParameters(2) ./ (xdata.^myParameters(2) + myParameters(3)) + myParameters(4));
fittedParameterslog=nlinfit(selectedxData,log(selectedyData),loghill,myParameters0)

%% plot etc

for plotindex = 1:2

    figure(plotindex), clf, hold on

    myXmin = min(allXData)/2;

    % plot fitted hill
    fitXData = logspace(log(myXmin)/log(10)/10,log(max(allXData))/log(10)*10,1000000);
    fittedHill = hill(fittedParameters,fitXData);
    plot(fitXData,fittedHill,'-','LineWidth',3,'Color',preferredcolors(2,:))
    % plot fitted loghill
    fittedHilllog = loghill(fittedParameterslog,fitXData);
    %figure
    plot(fitXData,exp(fittedHilllog),'-','LineWidth',3,'Color',preferredcolors(3,:))

    % plot data
    l=plot(allXData,allYData,'ok');
    set(l,'LineWidth',3);

    % plot no induction datapoint
    l=plot(myXmin,noinduction(2),'or');
    set(l,'LineWidth',3);

    % For user convenience make string that contains function in text
    hillEquationString = 'ymax * x^n / (x^n + K)';
    paramsString = ...
        ['{ymax=' num2str(fittedParameters(1),'%3.2s'),...
        ', n=' num2str(fittedParameters(1),'%.2e'),...
        ', K=' num2str(fittedParameters(2),'%.2e'),...
        ', offset=' num2str(fittedParameters(4),'%.2e'),...
        '}'];    

    % labels etc.
    xlabel('aTc (ng/ml)')
    ylabel('Luciferase activity (a.u.)')
    MW_makeplotlookbetter(16)
    myTitlepart1='ptet dose response curve (Lutz et al., 1997)'; 
    title([myTitlepart1 10 'Fitted (lin) Hill: ' hillEquationString, 10, paramsString],'Interpreter','None')
    title([myTitlepart1 10 'Fitted Hill (linear/log): ' hillEquationString, 10, 'scriptname: lutz1997plot.m'],'Interpreter','None')
    
    xlim([myXmin,max(allXData)*2])
    ylim([myYlimmin,myYlimmax])

    % logscale
    if plotindex==1
        set(gca,'xscale','log','yscale','log');
    else
        xlim([myXmin,max(allXData)*1.1])
        ylim([min([allYData])/2,max(allYData)*1.5])
    end
    
end

%% 
NUMBEROFSAMPLES = 6;
ATCLOW = 2; % lowest desired AtC inducer concentration
ATCHIGH = 20; % highest desired AtC inducer concentration

figure(1)
myCheckConcentrations = logspace(log(ATCHIGH)/log(10),log(ATCLOW)/log(10),NUMBEROFSAMPLES)
factor = myCheckConcentrations(1)/myCheckConcentrations(2)
for sampleLineIndex=1:NUMBEROFSAMPLES
    plot([myCheckConcentrations(sampleLineIndex),myCheckConcentrations(sampleLineIndex)],...
        [myYlimmin,myYlimmax],'--','Color',[.5 .5 .5],'LineWidth',2)
end

%% Extremely wide range
NUMBEROFSAMPLES = 9;
ATCLOW = 2; % lowest desired AtC inducer concentration
ATCHIGH = 40000; % highest desired AtC inducer concentration

figure(1)
myCheckConcentrations = logspace(log(ATCHIGH)/log(10),log(ATCLOW)/log(10),NUMBEROFSAMPLES)
mat2str(myCheckConcentrations)
factor = myCheckConcentrations(1)/myCheckConcentrations(2)
for sampleLineIndex=1:NUMBEROFSAMPLES
    plot([myCheckConcentrations(sampleLineIndex),myCheckConcentrations(sampleLineIndex)],...
        [myYlimmin,myYlimmax],'--','Color',[0 0 .7],'LineWidth',2)
end


%% arbitrary hill function
%{
figure(2)
loglog(allXData,hill([3000,5,4,1],allXData),'-','LineWidth',3,'Color',[.5,.5,.5])
xlim([min(allXData)/2,max(allXData)*2])
ylim([min(allYData)/2,max(allYData)*2])
%}