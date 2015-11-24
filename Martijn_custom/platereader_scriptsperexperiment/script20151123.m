PLOTDIR = 'U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\2015_11_23\summaryPlots\'

%% growth rate controls
FIELDNAMEY = 'muValues';

figure(1); clf; hold on

% collect data
values=[]; valuesstd=[];
for i=1:numel(controls)      
    
    % raw
    plot(ones(1,numel(controls(i).output.(FIELDNAMEY)))*i,[controls(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',2,'Color',[.6 .6 .6])
    % mean
    currentMean = mean([controls(i).output.(FIELDNAMEY)]);
    plot(i,currentMean,'ok','LineWidth',3)        
    % std        
    currentStd = std([controls(i).output.(FIELDNAMEY)]);
    plot(ones(1,2)*i,[currentMean-currentStd , currentMean+currentStd],'-','LineWidth',3,'Color','k')            
    
end

set(gca, 'XTickLabel',{controls.name}, 'XTick',1:numel(controls));

% cosmetics
theXlim = [0,numel(controls)+1];
xlim(theXlim)
ylim([0,1])
ylabel('Fitted growth rate [dbl/hr]')
MW_makeplotlookbetter(15);
plot(theXlim,[0,0],'k-')

%% fluor data controls
FIELDNAMEY = 'meanFluor';

figure(2); clf; hold on



% collect data
values=[]; valuesstd=[];
for i=1:numel(controls)      
    
    % raw
    plot(ones(1,numel(controls(i).output.(FIELDNAMEY)))*i,[controls(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',2,'Color',[.6 .6 .6])
    % mean
    currentMean = mean([controls(i).output.(FIELDNAMEY)]);
    plot(i,currentMean,'ok','LineWidth',3)        
    % std        
    currentStd = std([controls(i).output.(FIELDNAMEY)]);
    plot(ones(1,2)*i,[currentMean-currentStd , currentMean+currentStd],'-','LineWidth',3,'Color','k')            
    
end

set(gca, 'XTickLabel',{controls.name}, 'XTick',1:numel(controls));

% cosmetics
theXlim = [0,numel(controls)+1];
xlim(theXlim)
ylim([-100000,300000])
ylabel('Fluor/OD [a.u.]')
MW_makeplotlookbetter(15);
plot(theXlim,[0,0],'k-')

%% fluor vs. aTc (inducer)
FIELDNAMEY = 'muValues';
FIELDNAMEX = 'aTc';
MYXLIM=[10^0,40000];
MYYLIM =[0,1];
USEFORLEGEND=[1,3,5,6];

thedata=aTcoutput;

figure(1); clf; hold on

% collect data
values=[]; valuesstd=[]; lines=[];
for i=1:numel(thedata)      
    
    % raw
    l=plot([thedata(i).output.(FIELDNAMEX)],[thedata(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',3,'Color',[thedata(i).COLOR])    
    lines(end+1)=l;
    
    % since we use log-scale, also plot the 0 values at the left xlim
    myunplottableidxs = find([thedata(i).output.(FIELDNAMEX)]==0)
    plot(ones(1,numel(myunplottableidxs))*MYXLIM(1),[thedata(i).output(myunplottableidxs).(FIELDNAMEY)],'o','MarkerSize',20,'LineWidth',3,'Color',[thedata(i).COLOR])    
    
end

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
ylim(MYYLIM)
ylabel('Growth rate [dbl/hr]')
xlabel('Inducer aTc [ng/ml]')
MW_makeplotlookbetter(20);
plot(theXlim,[0,0],'k-')
title('script20151123.m');

legend(lines(USEFORLEGEND),{thedata(USEFORLEGEND).name},'Location','northeastoutside');

saveas(gcf,[PLOTDIR 'aTc_growth.png'],'png')
saveas(gcf,[PLOTDIR 'aTc_growth.eps'],'epsc')


%% fluor vs. aTc (inducer)
FIELDNAMEY = 'meanFluor';
FIELDNAMEX = 'aTc';
MYXLIM=[10^0,40000];
MYYLIM =[0,300000];
USEFORLEGEND=[1,3,5,6];

thedata=aTcoutput;

figure(2); clf; hold on

% collect data
values=[]; valuesstd=[]; lines=[];
for i=1:numel(thedata)      
    
    % raw
    l=plot([thedata(i).output.(FIELDNAMEX)],[thedata(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',3,'Color',[thedata(i).COLOR])    
    lines(end+1)=l;
    
    % since we use log-scale, also plot the 0 values at the left xlim
    myunplottableidxs = find([thedata(i).output.(FIELDNAMEX)]==0)
    plot(ones(1,numel(myunplottableidxs))*MYXLIM(1),[thedata(i).output(myunplottableidxs).(FIELDNAMEY)],'o','MarkerSize',20,'LineWidth',3,'Color',[thedata(i).COLOR])    
    
end

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
ylim(MYYLIM)
ylabel('Fluor reporter signal/OD [a.u.]')
xlabel('Inducer aTc [ng/ml]')
MW_makeplotlookbetter(20);
plot(theXlim,[0,0],'k-')
title('script20151123.m');

legend(lines(USEFORLEGEND),{thedata(USEFORLEGEND).name},'Location','northeastoutside');

saveas(gcf,[PLOTDIR 'aTc_fluor.png'],'png')
saveas(gcf,[PLOTDIR 'aTc_fluor.eps'],'epsc')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% growth rate cAMP extracell. experiment
FIELDNAMEY = 'muValues';

figure(1); clf; hold on

% collect data
values=[]; valuesstd=[];
for i=1:numel(alloutput)      
    
    % raw
    plot(ones(1,numel(alloutput(i).output.(FIELDNAMEY)))*i,[alloutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',2,'Color',[.6 .6 .6])
    % mean
    currentMean = mean([alloutput(i).output.(FIELDNAMEY)]);
    plot(i,currentMean,'ok','LineWidth',3)        
    % std        
    currentStd = std([alloutput(i).output.(FIELDNAMEY)]);
    plot(ones(1,2)*i,[currentMean-currentStd , currentMean+currentStd],'-','LineWidth',3,'Color','k')            
    
end

set(gca, 'XTickLabel',{alloutput.name}, 'XTick',1:numel(alloutput));

% cosmetics
theXlim = [0,numel(alloutput)+1];
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
for i=1:numel(alloutput)      
    
    % raw
    plot(ones(1,numel(alloutput(i).output.(FIELDNAMEY)))*i,[alloutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',2,'Color',[.6 .6 .6])
    % mean
    currentMean = mean([alloutput(i).output.(FIELDNAMEY)]);
    plot(i,currentMean,'ok','LineWidth',3)        
    % std        
    currentStd = std([alloutput(i).output.(FIELDNAMEY)]);
    plot(ones(1,2)*i,[currentMean-currentStd , currentMean+currentStd],'-','LineWidth',3,'Color','k')            
    
end

set(gca, 'XTickLabel',{alloutput.name}, 'XTick',1:numel(alloutput));

% cosmetics
theXlim = [0,numel(alloutput)+1];
xlim(theXlim)
ylim([-100000,300000])
ylabel('Fluor/OD [a.u.]')
MW_makeplotlookbetter(15);
plot(theXlim,[0,0],'k-')


%% growth rate vs. cAMP
FIELDNAMEY = 'muValues';
FIELDNAMEX = 'cAMP';
MYXLIM=[10^0,10000];
USEFORLEGEND=[1,3,5,6];

figure(1); clf; hold on

% collect data
values=[]; valuesstd=[]; lines=[];
for i=1:numel(alloutput)      
    
    % raw
    l=plot([alloutput(i).output.(FIELDNAMEX)],[alloutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',3,'Color',[alloutput(i).COLOR])    
    lines(end+1)=l;
    
    % since we use log-scale, also plot the 0 values at the left xlim
    myunplottableidxs = find([alloutput(i).output.(FIELDNAMEX)]==0)
    plot(ones(1,numel(myunplottableidxs))*MYXLIM(1),[alloutput(i).output(myunplottableidxs).(FIELDNAMEY)],'o','MarkerSize',25,'LineWidth',3,'Color',[alloutput(i).COLOR])    
    
end

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
ylim([0,1])
ylabel('Fitted growth rate [dbl/hr]')
xlabel('Extracellular cAMP [uM]')
MW_makeplotlookbetter(20);
plot(theXlim,[0,0],'k-')
title('script20151123.m');

legend(lines(USEFORLEGEND),{alloutput(USEFORLEGEND).name},'Location','northeastoutside');

saveas(gcf,[PLOTDIR 'cAMP_growth.png'],'png')
saveas(gcf,[PLOTDIR 'cAMP_growth.eps'],'epsc')

%% fluor vs. cAMP (extracellular cAMP)
FIELDNAMEY = 'meanFluor';
FIELDNAMEX = 'cAMP';
MYXLIM=[10^0,10000];
MYYLIM=[0 300000];
USEFORLEGEND=[1,3,5,6];

figure(1); clf; hold on

% collect data
values=[]; valuesstd=[]; lines=[];
for i=1:numel(alloutput)      
    
    % raw
    l=plot([alloutput(i).output.(FIELDNAMEX)],[alloutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',3,'Color',[alloutput(i).COLOR])    
    lines(end+1)=l;
    
    % since we use log-scale, also plot the 0 values at the left xlim
    myunplottableidxs = find([alloutput(i).output.(FIELDNAMEX)]==0)
    plot(ones(1,numel(myunplottableidxs))*MYXLIM(1),[alloutput(i).output(myunplottableidxs).(FIELDNAMEY)],'o','MarkerSize',25,'LineWidth',3,'Color',[alloutput(i).COLOR])    
    
end

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
ylim(MYYLIM)
ylabel('Fluor reporter signal/OD [a.u.]')
xlabel('Extracellular cAMP [uM]')
MW_makeplotlookbetter(20);
plot(theXlim,[0,0],'k-')
title('script20151123.m');

legend(lines(USEFORLEGEND),{alloutput(USEFORLEGEND).name},'Location','northeastoutside');

saveas(gcf,[PLOTDIR 'cAMP_fluor.png'],'png')
saveas(gcf,[PLOTDIR 'cAMP_fluor.eps'],'epsc')

%% scatter  fluor vs. growth (extracellular cAMP)
FIELDNAMEX = 'meanFluor';
FIELDNAMEY = 'muValues';
MYXLIM=[150,300000];
MYYLIM=[0 1];
USEFORLEGEND=[1,3,5,6];

figure(3); clf; hold on

% collect data
values=[]; valuesstd=[]; lines=[];
for i=[1:4]
    
    % raw
    l=plot([alloutput(i).output.(FIELDNAMEX)],[alloutput(i).output.(FIELDNAMEY)],'x','MarkerSize',25,'LineWidth',3,'Color',[alloutput(i).COLOR])    
    lines(end+1)=l;
    
    % since we use log-scale, also plot the 0 values at the left xlim
    myunplottableidxs = find([alloutput(i).output.(FIELDNAMEX)]==0)
    plot(ones(1,numel(myunplottableidxs))*MYXLIM(1),[alloutput(i).output(myunplottableidxs).(FIELDNAMEY)],'o','MarkerSize',25,'LineWidth',3,'Color',[alloutput(i).COLOR])    
    
end

% cosmetics
set(gca,'xscale','log');
xlim(MYXLIM)
ylim(MYYLIM)
ylabel('Growth rate [dbl/hr]')
xlabel('Fluor reporter signal/OD [a.u.]')
MW_makeplotlookbetter(20);
plot(theXlim,[0,0],'k-')
title('script20151123.m');

legend(lines(USEFORLEGEND),{alloutput(USEFORLEGEND).name},'Location','northeastoutside');

saveas(gcf,[PLOTDIR 'cAMP_scatter.png'],'png')
saveas(gcf,[PLOTDIR 'cAMP_scatter.eps'],'epsc')

