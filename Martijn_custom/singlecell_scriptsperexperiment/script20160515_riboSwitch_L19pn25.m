
% Script to plot effect of ribosomal switch.
% Script adapted from script20160126_plottingSarahsData.
%
%
%
%
% PS. See general_playingwithcorrelations.m for some notes on coefficients.

SCRIPTNAME = 'script20160515_riboSwitch_L19pn25.m';

%% Loading datasets

%% pn25-mCherry
% \\storage01\data\AMOLF\groups\tans-group\Former-Users\Boulineau\Data analysis\Antibiotics\2012-07-12
% shift @ 366 min (pos 315), .5uM
TITLE='pn25-mCherry';
THESHIFTTIME=366;
EXPERIMENTALID = ['2012-07-12_pos6_' TITLE];
GROWTHFIELD = 'muP21_fitNew_all';
load('\\storage01\data\AMOLF\groups\tans-group\Former-Users\Boulineau\Data analysis\Antibiotics\2012-07-12\pos6crop\data\pos6crop-Schnitz.mat');
subdir = [EXPERIMENTALID '\']

%% 1 uM
TITLE='L31-mCherry';
THESHIFTTIME=309;
EXPERIMENTALID = ['2012-11-15_pos5_' TITLE];
GROWTHFIELD = 'muP17_fitNew_all';
load('\\storage01\data\AMOLF\groups\tans-group\Former-Users\Boulineau\Data analysis\Antibiotics\2012-11-15\pos5crop\data\pos5crop-Schnitz.mat');
subdir = [EXPERIMENTALID '\']

%% .29 uM
TITLE = 'L31-mCerulean'
THESHIFTTIME = 171.1833; 
% timeofswitch = datenum(2016,02,11,01,20,00); schnitzstarttime = min([schnitzcells.timestamp]); switchedafterdays=timeofswitch-schnitzstarttime; switchedafterminutes=switchedafterdays*24*60
EXPERIMENTALID = ['2016-02-10_pos5_' TITLE];
GROWTHFIELD = 'muP5_fitNew_all';
load('F:\A_Tans1_step1_incoming_not_backed_up\2016-02-10\pos5\data\pos5-Schnitz.mat');
subdir = [EXPERIMENTALID '\']

%% .

TITLE = 'L31-mCerulean'
THESHIFTTIME = 226.9833; % Note 221 mins is what you'd manually calculate.
% timeofswitch = datenum(2016,02,17,16,23,00); schnitzstarttime = min([schnitzcells.timestamp]); switchedafterdays=timeofswitch-schnitzstarttime; switchedafterminutes=switchedafterdays*24*60
EXPERIMENTALID = ['2016-02-17_pos2_' TITLE];
GROWTHFIELD = 'muP5_fitNew_all';
load('F:\A_Tans1_step1_incoming_not_backed_up\2016-02-17\pos2\data\pos2-Schnitz.mat');
subdir = [EXPERIMENTALID '\']

%% 

TITLE = 'L19-mCerulean_pn25-yfp'
THESHIFTTIME = (24*60+46)-(22*60+9); % Note 221 mins is what you'd manually calculate.
% timeofswitch = datenum(2016,02,17,16,23,00); schnitzstarttime = min([schnitzcells.timestamp]); switchedafterdays=timeofswitch-schnitzstarttime; switchedafterminutes=switchedafterdays*24*60
EXPERIMENTALID = ['2016-05-15_pos1_' TITLE];
GROWTHFIELD = 'muP5_fitNew_all';
load('F:\A_Tans1_step1_incoming_not_backed_up\2016-02-17\pos2\data\pos2-Schnitz.mat');
subdir = [EXPERIMENTALID '\']


%% If necessary
%schnitzcells = s_rm;

%% Parameters
OUTPUTDIR = 'D:\Local_Software\Martijn_extensions\Martijn_custom\singlecell_scriptsperexperiment\plots20160217\';

if ~exist(OUTPUTDIR,'dir')
    mkdir(OUTPUTDIR)
end

if ~exist([OUTPUTDIR subdir],'dir')
    mkdir([OUTPUTDIR subdir])
end

finalTitle = [TITLE 10 SCRIPTNAME];

%% backwards compatibility
if ~isfield(schnitzcells,'frame_nrs')
    schnitzcells = MW_calculateframe_nrs(schnitzcells);
end

% paramters
allFrameNrs = unique([schnitzcells.frame_nrs]);
minFrameNr = min(allFrameNrs);
maxFrameNr = max(allFrameNrs);

% obtain growth per frame number
timeField = 'time';
fieldOfInterest = GROWTHFIELD;
[growthtimesPerframe_nrs, growthPerframe_nrs, growthschnitzesPerFrame_nrs, growthplottingfXaxis] = ...
    MW_getparamfromschnitzcells(schnitzcells, timeField, fieldOfInterest)

% obtain fluor per frame number
timeField = 'time';
fieldOfInterest = 'C5_mean_all';
[mcherryTimesPerframe_nrs, mcherryPerframe_nrs, mcherryschnitzesPerFrame_nrs, mcherryplottingfXaxis] = ...
    MW_getparamfromschnitzcells(schnitzcells, timeField, fieldOfInterest)


%% Plot some of the data

%% Some additional parameters

somecolors = linspecer(3);

% Determine frames with data (brute force)
IdxWithData = [];
for i = 1:numel(mcherryPerframe_nrs)
    if ~any(isnan(mcherryPerframe_nrs{i}))
        IdxWithData(end+1) = i;
    end
end
maxIdxWithData = max(IdxWithData);

% Get # datapoints per frame
% Get time at which frame was taken
numelgrowthtimesPerframe_nrs = numel(growthtimesPerframe_nrs);
dataCount = []; allTimes = NaN(1,numelgrowthtimesPerframe_nrs);
for i = 1:numel(growthtimesPerframe_nrs)
    dataCount(i) = numel(growthtimesPerframe_nrs{i});
    
    if any(~(growthtimesPerframe_nrs{i}==growthtimesPerframe_nrs{i}(1)))
        error('Unequal timefields! Assumptions broken.');
        % The output for the timefield correspond to each schnitz that is a
        % member of this frame, these times should all be identical as they
        % come from the same frame. Hence I can use the first entry. If
        % note, something went horribly wrong ;).
    end
    allTimes(i) = growthtimesPerframe_nrs{i}(1);
end

minTime = min(allTimes);
maxTime = max(allTimes);

%% average growth
figure; clf; hold on;
meangrowth = [];
for i = 1:numel(growthtimesPerframe_nrs)
    meangrowth(i) = mean(growthPerframe_nrs{i});
    stdgrowth(i) = std(growthPerframe_nrs{i});
    
    %plot(growthtimesPerframe_nrs{i}(1),meangrowth(i),'.'); 
    plot(growthtimesPerframe_nrs{i},growthPerframe_nrs{i},'xb','Color',[.5 .5 .5]); 
    plot(growthtimesPerframe_nrs{i}(1),meangrowth(i),'.k'); 
    %errorbar(growthtimesPerframe_nrs{i}(1),meangrowth(i),stdgrowth(i),'k.'); 
end

%% single cell growth 
figure; clf; hold on;
for i = 1:numel(growthtimesPerframe_nrs)
    plot(growthtimesPerframe_nrs{i},growthPerframe_nrs{i},'.')
end

%% total mcherry
figure; clf; hold on;
for i = 1:numel(mcherryTimesPerframe_nrs)
    plot(mcherryTimesPerframe_nrs{i}(1),sum(mcherryPerframe_nrs{i}),'.r')
end

%% single cell mcherry 
figure; clf; hold on;
for i = 1:numel(mcherryTimesPerframe_nrs)
    plot(mcherryTimesPerframe_nrs{i},mcherryPerframe_nrs{i},'.r')
end

%% scatter for each frame, plus extract correlation
%frameNr = 502;

betaYs=NaN(1,numelgrowthtimesPerframe_nrs);
fitcoefficientsa=NaN(1,numelgrowthtimesPerframe_nrs);
fitcoefficientsb=NaN(1,numelgrowthtimesPerframe_nrs);
RXY=NaN(1,numelgrowthtimesPerframe_nrs);
for Idx = IdxWithData % note frameIdx ~= frameNr
    % scatter figure
    h=figure(100); clf; hold on;
    plot(mcherryPerframe_nrs{Idx},growthPerframe_nrs{Idx},'.')
    
    % calculate regression coefficient
    C = cov(mcherryPerframe_nrs{Idx},growthPerframe_nrs{Idx});
    betaYs(Idx) = C(1,2)/var(mcherryPerframe_nrs{Idx});
    RXY(Idx) = C(1,2)/sqrt(var(mcherryPerframe_nrs{Idx})*var(growthPerframe_nrs{Idx}))
    
    % calculate linear fit
    [p,S] = polyfit(mcherryPerframe_nrs{Idx},growthPerframe_nrs{Idx},1);
    fitcoefficientsa(Idx) = p(2);
    fitcoefficientsb(Idx) = p(1);
    % plot fitted slope
    fitxlim = [min(mcherryPerframe_nrs{Idx}), max(mcherryPerframe_nrs{Idx})];
    fitx = fitxlim;  % need only two x points since line
    fity = fitxlim*p(1)+p(2);
    plot(fitx, fity);
        % note that an indication of reliability is lacking.
    
    xlabel('mCherry signal [a.u.]');
    ylabel('Growth rate [dbl/hr]');
        
    title(['Frame ' num2str(allFrameNrs(Idx)) ', time = ' num2str(mcherryTimesPerframe_nrs{Idx}(1)) ...
            ', slope = ' num2str(betaYs(Idx)) ', cells = ' num2str(dataCount(Idx)),...
            10 'correlation R(growth, mCherry) = ' num2str(RXY(Idx))]);
    saveas(h, [OUTPUTDIR subdir 'scatter_fr' sprintf('%04d',Idx) '.tif'],'tif');
end

%% make plot with overview subplots

h=figure(1); clf;
MYXLIM=[minTime,maxTime];
MYYLIM1 = [1, max(dataCount)];%*10];
%MYYLIM1 = [1, max(dataCount)*1.1];
MYYLIM2 = [0, max(meangrowth)];
MYYLIM3 = [min(betaYs), max(betaYs)];

% subplot 1
subplot(4,1,1); hold on;
plot([THESHIFTTIME, THESHIFTTIME],MYYLIM1,'k-','LineWidth',2); % shift time
l=plot(allTimes,dataCount); % actual data
%[a,l,l2]=plotyy(allTimes,dataCount,allTimes,meangrowth); % data
% log scale
%set(gca,'yscale','log'); 
xlim(MYXLIM); ylim(MYYLIM1);
% cosmetics
MW_makeplotlookbetter(15);
set(l, 'LineWidth', 3, 'Color', 'k');%somecolors(1,:));
ylabel('# cells')

title(finalTitle,'Interpreter','None');

% subplot 2 - growth rates
subplot(4,1,2); hold on;
%== data
%l=plot(allTimes,meangrowth); % data
for i = 1:numel(growthtimesPerframe_nrs)
    plot(growthtimesPerframe_nrs{i},growthPerframe_nrs{i},'xb','Color',[.5 .5 .5]);     
end
l = plot(allTimes,meangrowth,'-k');
%==
xlim(MYXLIM); ylim(MYYLIM2);
% cosmetics
MW_makeplotlookbetter(15);
set(l, 'LineWidth', 3, 'Color', 'k');%somecolors(1,:));
ylabel('dbl/hr')
% shift time
plot([THESHIFTTIME, THESHIFTTIME],MYYLIM2,'k-','LineWidth',2); % shift time


% subplot 3 (correlation coefficient)
subplot(4,1,3); hold on;
plot(MYXLIM,[0,0],'k-')% 0 axis
plot([THESHIFTTIME, THESHIFTTIME],MYYLIM3,'k-','LineWidth',2); % shift time
l=plot(allTimes, fitcoefficientsb,'x'); % fit coefficient
xlim(MYXLIM); ylim(MYYLIM3);
% cosmetics
MW_makeplotlookbetter(15);
set(l, 'LineWidth', 3, 'Color', 'k');%somecolors(1,:));
xlabel('Time [min]');
ylabel('Slope');


% subplot 4 (correlation coefficient)
subplot(4,1,4); hold on;
plot(MYXLIM,[0,0],'k-')% 0 axis
plot([THESHIFTTIME, THESHIFTTIME],[-1,1],'k-','LineWidth',2); % shift time
l=plot(allTimes, RXY,'x'); % fit coefficient
xlim(MYXLIM); ylim([-1,1]);
% cosmetics
MW_makeplotlookbetter(15);
set(l, 'LineWidth', 3, 'Color', 'k');%somecolors(1,:));
xlabel('Time [min]');
ylabel('R(ribo,growth)');

saveas(h, [OUTPUTDIR subdir 'summaryPlot.tif']);
saveas(h, [OUTPUTDIR subdir 'summaryPlot.fig']);
saveas(h, [OUTPUTDIR subdir 'summaryPlot.eps'],'epsc');


%% Now plot clouds

% median-normalized

figure(2); clf;

numelIdxWithData = numel(IdxWithData);
sizeSubplot = ceil(sqrt(numelIdxWithData+extraplots));

% axis label plot
extraplots = 1;
subplottightmargin(sizeSubplot,sizeSubplot,1,0.04); hold on;
set(gca,'xtick',[],'ytick',[])
plot([-1,1],[0,0],'k-');
plot([0,0],[-1,1],'k-');
text(.5,-0.5,'\rightarrow ribo')
text(.1,.7,'\uparrow growth')

xlim([-1.5,1.5]);
ylim([-2,1.5]);
        
% go over subplots
theColors = linspecer(2);
for y = 1:sizeSubplot
    for x = 1:sizeSubplot   
        
        % calculate which subplot
        subplotIdx = (y-1)*sizeSubplot+x; 
            
        % Get idx
        Idx = IdxWithData(subplotIdx);        
        
        % stop plotting if nothing to plot
        if subplotIdx > numelIdxWithData; break; end % ugly
        
        % determine time & color of plot
        % ===
        currentTime = mcherryTimesPerframe_nrs{Idx}(1);        
        % color
        if THESHIFTTIME<currentTime            
            currentColor = theColors(1,:);
        else
            currentColor = theColors(2,:);
        end
        
        % create plot
        %subplot(sizeSubplot,sizeSubplot,subplotIdx); hold on;
        subplottightmargin(sizeSubplot,sizeSubplot,subplotIdx+extraplots,0.04); hold on;
        set(gca,'xtick',[],'ytick',[])

        % plot cloud        
        xdata = mcherryPerframe_nrs{Idx};
        ydata = growthPerframe_nrs{Idx};
        plot(   (xdata-median(xdata))/median(xdata),...
                (ydata-median(ydata))/median(ydata),...,
                '.', 'Color',currentColor)
        plot([-1,1],[0,0],'k-');
        plot([0,0],[-1,1],'k-');
        
        xlim([-1.5,1.5]);
        ylim([-2,1.5]);
        
        plotinfo = ['R=' sprintf('%0.2f', RXY(Idx)) ', t=' sprintf('%3.0f', allTimes(Idx))];
        
        text(0,-1.75,plotinfo,'HorizontalAlignment','Center')
                
    end
end

%set(gcf,'Units','normal')
%set(gca,'Position',[0 0 1 1])













