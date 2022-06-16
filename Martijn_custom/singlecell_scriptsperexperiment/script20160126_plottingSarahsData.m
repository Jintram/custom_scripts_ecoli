

% Probably the usefull datasets are the .5uM tetracycline switches. There
% should be datasets that have the ASC666 strain, which has L31-mCherry
% ribosomal label.
%
% Some of Sarah's data seems to be preliminary analysis only.
%
%
% PS. See general_playingwithcorrelations.m for some notes on coefficients.

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


%% If necessary
%schnitzcells = s_rm;

%% Parameters
OUTPUTDIR = 'D:\Local_Software\Martijn_extensions\Martijn_custom\singlecell_scriptsperexperiment\plots20160126\';

if ~exist(OUTPUTDIR,'dir')
    mkdir(OUTPUTDIR)
end

if ~exist([OUTPUTDIR subdir],'dir')
    mkdir([OUTPUTDIR subdir])
end


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
dataCount = []; allTimes = NaN(1,numel(growthtimesPerframe_nrs));
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
    
    plot(growthtimesPerframe_nrs{i}(1),meangrowth(i),'.'); 
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

betaYs=NaN(1,maxIdxWithData);
fitcoefficientsa=NaN(1,maxIdxWithData);
fitcoefficientsb=NaN(1,maxIdxWithData);
for Idx = IdxWithData % note frameIdx ~= frameNr
    % scatter figure
    h=figure(100); clf; hold on;
    plot(mcherryPerframe_nrs{Idx},growthPerframe_nrs{Idx},'.')
    
    % calculate regression coefficient
    C = cov(mcherryPerframe_nrs{Idx},growthPerframe_nrs{Idx});
    betaYs(Idx) = C(1,2)/var(mcherryPerframe_nrs{Idx});
    
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
    
    title(['Frame ' num2str(allFrameNrs(Idx)) ', time = ' num2str(mcherryTimesPerframe_nrs{Idx}(1)) ...
            ', slope = ' num2str(betaYs(Idx)) ', cells = ' num2str(dataCount(Idx))]);
    saveas(h, [OUTPUTDIR subdir 'scatter_fr' sprintf('%04d',Idx) '.tif'],'tif');
end

%% 
h=figure; clf;
MYXLIM=[minTime,maxTime];
MYYLIM1 = [1, max(dataCount)*10];
%MYYLIM1 = [1, max(dataCount)*1.1];
MYYLIM2 = [0, max(meangrowth)];
MYYLIM3 = [min(betaYs), max(betaYs)];

% subplot 1
subplot(3,1,1); hold on;
plot([THESHIFTTIME, THESHIFTTIME],MYYLIM1,'k:'); % shift time
l=plot(allTimes,dataCount); % actual data
%[a,l,l2]=plotyy(allTimes,dataCount,allTimes,meangrowth); % data
set(gca,'yscale','log');
xlim(MYXLIM); ylim(MYYLIM1);
% cosmetics
MW_makeplotlookbetter(15);
set(l, 'LineWidth', 3, 'Color', somecolors(1,:));
ylabel('# cells')
title(TITLE);


% subplot 2
subplot(3,1,2); hold on;
plot([THESHIFTTIME, THESHIFTTIME],MYYLIM2,'k:'); % shift time
l=plot(allTimes,meangrowth); % data
xlim(MYXLIM); ylim(MYYLIM2);
% cosmetics
MW_makeplotlookbetter(15);
set(l, 'LineWidth', 3, 'Color', somecolors(1,:));
ylabel('dbl/hr')


% subplot 3
subplot(3,1,3); hold on;
plot(MYXLIM,[0,0],'k-')% 0 axis
plot([THESHIFTTIME, THESHIFTTIME],MYYLIM3,'k:'); % shift time
l=plot(allTimes, fitcoefficientsb,'x'); % fit coefficient
xlim(MYXLIM); ylim(MYYLIM3);
% cosmetics
MW_makeplotlookbetter(15);
set(l, 'LineWidth', 3, 'Color', somecolors(1,:));
xlabel('Time [min]');
ylabel('Regression coef.');

saveas(h, [OUTPUTDIR subdir 'summaryPlot.tif']);
saveas(h, [OUTPUTDIR subdir 'summaryPlot.fig']);
saveas(h, [OUTPUTDIR subdir 'summaryPlot.eps'],'epsc');



