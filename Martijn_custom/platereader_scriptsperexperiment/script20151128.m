
%% Load the dataset

load('U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\2015_11_28\alloutput.mat');
PLOTDIR = 'U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\2015_11_28\summaryPlots\rerun\'

%% Converting names to cAMP concentration values.

for i = 1:numel(ampoutput)
    for j = 1:numel(ampoutput(i).output)
        ampoutput(i).output(j).cAMP = str2double(ampoutput(i).output(j).groupName(4:8))
    end
end

ampoutput(1).COLOR=[1 0 0];
ampoutput(2).COLOR=[1 0 0];
ampoutput(3).COLOR=[0 0 1];
ampoutput(4).COLOR=[0 0 1];
    
%% growth rate cAMP extracell. experiment
FIELDNAMEY = 'muValues';

figure(1); clf; hold on

% collect data
values=[]; valuesstd=[];
for i=1:numel(ampoutput)      
    
    % raw
    plot(ones(1,numel(ampoutput(i).output.(FIELDNAMEY)))*i,[ampoutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',2,'Color',[.6 .6 .6])
    % mean
    currentMean = mean([ampoutput(i).output.(FIELDNAMEY)]);
    plot(i,currentMean,'ok','LineWidth',3)        
    % std        
    currentStd = std([ampoutput(i).output.(FIELDNAMEY)]);
    plot(ones(1,2)*i,[currentMean-currentStd , currentMean+currentStd],'-','LineWidth',3,'Color','k')            
    
end

set(gca, 'XTickLabel',{ampoutput.name}, 'XTick',1:numel(ampoutput));

% cosmetics
theXlim = [0,numel(ampoutput)+1];
xlim(theXlim)
ylim([0,1])
ylabel('Fitted growth rate [dbl/hr]')
MW_makeplotlookbetter(15);
plot(theXlim,[0,0],'k-')

%% fluor data output
FIELDNAMEY = 'meanFluor';

figure(2); clf; hold on

% collect data
values=[]; valuesstd=[];
for i=1:numel(ampoutput)      
    
    % raw
    plot(ones(1,numel(ampoutput(i).output.(FIELDNAMEY)))*i,[ampoutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',2,'Color',[.6 .6 .6])
    % mean
    currentMean = mean([ampoutput(i).output.(FIELDNAMEY)]);
    plot(i,currentMean,'ok','LineWidth',3)        
    % std        
    currentStd = std([ampoutput(i).output.(FIELDNAMEY)]);
    plot(ones(1,2)*i,[currentMean-currentStd , currentMean+currentStd],'-','LineWidth',3,'Color','k')            
    
end

set(gca, 'XTickLabel',{ampoutput.name}, 'XTick',1:numel(ampoutput));

% cosmetics
theXlim = [0,numel(ampoutput)+1];
xlim(theXlim)
ylim([-100000,300000])
ylabel('Fluor/OD [a.u.]')
MW_makeplotlookbetter(15);
plot(theXlim,[0,0],'k-')


%% growth rate vs. cAMP
FIELDNAMEY = 'muValues';
FIELDNAMEX = 'cAMP';
MYXLIM=[10^0,1000000];
USEFORLEGEND=[1,3];

% For later use
if ~exist('xydatacollection','var'), xydatacollection = []; end

figure(1); 
if ~exist('NOCLEARFIGURE','var'), clf; end;
hold on

% collect data
values=[]; valuesstd=[]; lines=[];
for i=1:numel(ampoutput)      
    
    % raw
    l=plot([ampoutput(i).output.(FIELDNAMEX)],[ampoutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',3,'Color',[ampoutput(i).COLOR])    
    lines(end+1)=l;
    
    % since we use log-scale, also plot the 0 values at the left xlim
    myunplottableidxs = find([ampoutput(i).output.(FIELDNAMEX)]==0);
    plot(ones(1,numel(myunplottableidxs))*MYXLIM(1),[ampoutput(i).output(myunplottableidxs).(FIELDNAMEY)],'o','MarkerSize',25,'LineWidth',3,'Color',[ampoutput(i).COLOR])    
    
    % For later use outside script
    xydatacollection = [xydatacollection [[ampoutput(i).output.(FIELDNAMEX)];[ampoutput(i).output.(FIELDNAMEY)]]];
    
end

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
ylim([0,1])
ylabel('Fitted growth rate [dbl/hr]')
xlabel('Extracellular cAMP [uM]')
MW_makeplotlookbetter(20);
plot(theXlim,[0,0],'k-')
title('script20151128.m');

legend(lines(USEFORLEGEND),{ampoutput(USEFORLEGEND).name},'Location','northeastoutside');

saveas(gcf,[PLOTDIR 'cAMP_growth.png'],'png')
saveas(gcf,[PLOTDIR 'cAMP_growth.eps'],'epsc')

%% fluor vs. cAMP (extracellular cAMP)
FIELDNAMEY = 'meanFluor';
FIELDNAMEX = 'cAMP';
MYXLIM=[10^0,100000];
MYYLIM=[0 300000];
USEFORLEGEND=[1,3];

figure(1); clf; hold on

% collect data
values=[]; valuesstd=[]; lines=[];
for i=1:numel(ampoutput)      
    
    % raw
    l=plot([ampoutput(i).output.(FIELDNAMEX)],[ampoutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',3,'Color',[ampoutput(i).COLOR])    
    lines(end+1)=l;
    
    % since we use log-scale, also plot the 0 values at the left xlim
    myunplottableidxs = find([ampoutput(i).output.(FIELDNAMEX)]==0)
    plot(ones(1,numel(myunplottableidxs))*MYXLIM(1),[ampoutput(i).output(myunplottableidxs).(FIELDNAMEY)],'o','MarkerSize',25,'LineWidth',3,'Color',[ampoutput(i).COLOR])    
    
end

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
ylim(MYYLIM)
ylabel('Fluor reporter signal/OD [a.u.]')
xlabel('Extracellular cAMP [uM]')
MW_makeplotlookbetter(20);
plot(theXlim,[0,0],'k-')
title('script20151128.m');

legend(lines(USEFORLEGEND),{ampoutput(USEFORLEGEND).name},'Location','northeastoutside');

saveas(gcf,[PLOTDIR 'cAMP_fluor.png'],'png')
saveas(gcf,[PLOTDIR 'cAMP_fluor.eps'],'epsc')

%% scatter  fluor vs. growth (extracellular cAMP)
FIELDNAMEX = 'meanFluor';
FIELDNAMEY = 'muValues';
MYXLIM=[150,300000];
MYYLIM=[0 1];
USEFORLEGEND=[1,3];

figure(3); clf; hold on

% collect data
values=[]; valuesstd=[]; lines=[];
for i=[1:4]
    
    % raw
    l=plot([ampoutput(i).output.(FIELDNAMEX)],[ampoutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',3,'Color',[ampoutput(i).COLOR])    
    lines(end+1)=l;
    
    % since we use log-scale, also plot the 0 values at the left xlim
    myunplottableidxs = find([ampoutput(i).output.(FIELDNAMEX)]==0)
    plot(ones(1,numel(myunplottableidxs))*MYXLIM(1),[ampoutput(i).output(myunplottableidxs).(FIELDNAMEY)],'o','MarkerSize',25,'LineWidth',3,'Color',[ampoutput(i).COLOR])    
    
end

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
ylim(MYYLIM)
ylabel('Growth rate [dbl/hr]')
xlabel('Fluor reporter signal/OD [a.u.]')
MW_makeplotlookbetter(20);
plot(theXlim,[0,0],'k-')
title('script20151128.m');

legend(lines(USEFORLEGEND),{ampoutput(USEFORLEGEND).name},'Location','northeastoutside');

saveas(gcf,[PLOTDIR 'cAMP_scatter.png'],'png')
saveas(gcf,[PLOTDIR 'cAMP_scatter.eps'],'epsc')


%% Extra plot

% xydatacollection obtained from executing script20151128 and 
% script20151128 consecutively

MYXLIM = [10^0, 10^5];

figure(4); clf; hold on;
plot(xydatacollection(1,:),xydatacollection(2,:),'ok','LineWidth',3)

meanvaluearray = [];
for xvalue = unique(xydatacollection(1,:))
    indicesthisxvalue = find(xydatacollection(1,:)==xvalue);
    plot(xydatacollection(1,indicesthisxvalue),xydatacollection(2,indicesthisxvalue),'o','LineWidth',3,'Color',[.6 .6 .6]);
    if xvalue == 0        
        plot(ones(1,numel(xydatacollection(2,indicesthisxvalue))).*MYXLIM(1),xydatacollection(2,indicesthisxvalue),'o','LineWidth',3,'Color',[.6 .6 .6]);
    end
    meanvaluearray = [meanvaluearray, [mean(xydatacollection(1,indicesthisxvalue)); mean(xydatacollection(2,indicesthisxvalue))]];
end

plot(meanvaluearray(1,:), meanvaluearray(2,:),'o--','LineWidth',3,'Color','k')

zerovaluexindex = find(meanvaluearray(1,:)==0)
plot(MYXLIM(1), meanvaluearray(2,zerovaluexindex),'>--','LineWidth',3,'Color','k')

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
ylim([0,1])
ylabel('Fitted growth rate [dbl/hr]')
xlabel('Extracellular cAMP [uM]')
MW_makeplotlookbetter(20);

title('script20151128.m')

saveas(gcf,[PLOTDIR 'twoexperimentscombined-Ocurve.png'],'png')
saveas(gcf,[PLOTDIR 'twoexperimentscombined-Ocurve.eps'],'epsc')

%% Add microscope data to figure

figure(4); hold on;

microscopedata = ...
[ 314, .567;...
 314, .529;...
 314, .611;...
 314, .705;...
 314, .510;...
 314, .675;...
 314, .623;...
  50, .386;...
1300, .821;...
1300, .80;...
1300, .727;...
1300, .763;...
1300, .818;...
1300, .810;...
1300, .774;]

plot(microscopedata(:,1),microscopedata(:,2),'x','LineWidth',3,'Color','r','MarkerSize',15)

disp('done');

%% 

%% fitTime(1) (proxy for lag phase) vs. cAMP
FIELDNAMEX = 'cAMP';
FIELDNAMEY = 'fitTimes';
MYXLIM=[10^0,100000];
USEFORLEGEND=[1,3];

% For later use
if ~exist('xydatacollection','var'), xydatacollection = []; end

figure(5); clf; hold on;

% go over data
values=[]; valuesstd=[]; lines=[];
for i=1:numel(ampoutput)      
    
    % raw
    for j = 1:numel([ampoutput(i).output])
    
        if ampoutput(i).output(j).(FIELDNAMEX) == 0
            ydata=[ampoutput(i).output(j).(FIELDNAMEY)(1)];
            plot(MYXLIM(1),ydata,'o','MarkerSize',25,'LineWidth',3,'Color',[ampoutput(i).COLOR]);
        else        
            ydata=[ampoutput(i).output(j).(FIELDNAMEY)(1)];
            l=plot([ampoutput(i).output(j).(FIELDNAMEX)],ydata,'x','MarkerSize',25,'LineWidth',3,'Color',[ampoutput(i).COLOR]);            
        end
    end
    
    lines(end+1)=l;
        
end

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
%ylim([0,1])
ylabel('Treshold passage (hrs)')
xlabel('Extracellular cAMP [uM]')
MW_makeplotlookbetter(20);
plot(theXlim,[0,0],'k-')
title('script20151128.m');

legend(lines(USEFORLEGEND),{ampoutput(USEFORLEGEND).name},'Location','northeastoutside');

saveas(gcf,[PLOTDIR 'treshold.png'],'png')
saveas(gcf,[PLOTDIR 'treshold.eps'],'epsc')