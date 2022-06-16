
%% 2015/02 MW 
% 
% Looking again at some data Aileen took.
% Question is whether cells are in exponential phase.

% Paramater settings ----------
FITTIME = 100;
EXPORTFOLDER = 'U:\PROJECTS\D_Projects_Aileen\images\';
% -----------------------------

% List of files analyzed in 2008
myFile = 'F:\X_Other_datasets\2008mondate\2008-02-23\Pos09-mini-01\data\Pos09-mini-01-Schnitz.mat';
myFile = 'F:\X_Other_datasets\2008mondate\2008-02-23\Pos10-mini-01\data\Pos10-mini-01-Schnitz.mat';
myFile = 'F:\X_Other_datasets\2008mondate\2008-02-23\ONcells-mini-01\data\ONcells-mini-01-Schnitz.mat';
myFile = 'F:\X_Other_datasets\2008mondate\2008-02-13\starting-mini-01\data\starting-mini-01-Schnitz.mat';

% List of files analyzed in 2007
% ===
% 2007-01-23
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-01-23\Pos05-mini-01\data\Pos05-mini-01-Schnitz.mat';
% 2007-06-19
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-06-19\jumping cy Pos10-mini-01\data\Pos10-mini-01-Schnitz.mat';
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-06-19\Pos01-mini-01\data\Pos01-mini-01-Schnitz.mat';
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-06-19\Pos02-mini-01\data\Pos02-mini-01-Schnitz.mat';
% 2007-08-16
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-08-16\Pos06-mini-01\data\Pos06-mini-01-Schnitz.mat';
% 2007-08-23
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-08-23\Pos01-mini-01\data\Pos01-mini-01-Schnitz.mat';
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-08-23\Pos03-mini-01\data\Pos03-mini-01-Schnitz.mat';
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-08-23\Pos08-mini-01\data\Pos08-mini-01-Schnitz.mat';
% redundant for 2007-08-23
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-08-23\Pos03-mini-01\Pos03-mini-01-Schnitz.mat';
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-08-23\short w cher Pos08-mini-01\Pos08-mini-01-Schnitz.mat';
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-08-23\early Pos08-mini-01\data\Pos08-mini-01-Schnitz.mat';

%%% CURRENT ONE: %%%
myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-08-23\early Pos08-mini-01\data\Pos08-mini-01-Schnitz.mat';
%%% //////////// %%%

% create identifier
%slashes = find(myFile=='\');
%identifier = myFile(slashes(end-4):end);
identifier = regexprep(myFile,'\\','-');
identifier = regexprep(identifier,'\.','-');

load(myFile); 
myFile % just print for user convenience

timeField = 'mins';
lengthField = 'lengthMicrons';

%% Analyze

myFrames = unique([schnitzcells.frames]);
myTimes = unique([schnitzcells.(timeField)]);
allLengths = [schnitzcells.(lengthField)];

% Surely this can be done prettier, but let's just loop to get the desired
% data now.
lengthsPerFrames = {}; plottingfXaxis = {}; timesPerFrames = {};
for f = myFrames
   
    lengthForThisFrame = [];
    timesForThisFrame = [];
    
    for schnitzIdx = 1:numel(schnitzcells)
        % Look whether this schnitz lives in current frame number,
        % and what timepoint belongs to the current frame.
        pointInSchnitz = find(schnitzcells(schnitzIdx).frames==f);
        
        % If so, add length at that timepoint to the collection of lengths
        if ~isempty(pointInSchnitz)
            timesForThisFrame(end+1) = schnitzcells(schnitzIdx).(timeField)(pointInSchnitz);            
            lengthForThisFrame(end+1) = schnitzcells(schnitzIdx).(lengthField)(pointInSchnitz);
        end
    end
    
    timesPerFrames{end+1} = timesForThisFrame;
    lengthsPerFrames{end+1} = lengthForThisFrame;
    plottingfXaxis{end+1} = ones(numel(lengthForThisFrame),1)*f; % convenient as x-axis later
    
end

%% Plotting

% Some plotting settings
h = figure('Name', myFile,'pos', [100 100 800+100 400+100]); % 700x600
subplot(1,2,1);
title([ 'Raw data, lengths bacteria at each frame'])
xlabel('Time in minutes');
ylabel('Length individual bacteria ({\mu}m)');
hold on;
xlim([min(myTimes),max(myTimes)])
ylim([min(allLengths),max(allLengths)])

summedLengthsPerFrame=[];
for i = 1:numel(plottingfXaxis)
    
    % raw data
    xdata = cell2mat(timesPerFrames(i));
    ydata = cell2mat(lengthsPerFrames(i));    
    plot(xdata, ydata,'o')
    set(gca,'FontSize',15)
        
    % determine sum of cell lengths
    summedLengthsPerFrame(end+1) = sum(cell2mat(lengthsPerFrames(i)));
    
end

% Plotting setup
subplot(1,2,2); 
semilogy(myTimes, summedLengthsPerFrame,'o');
hold on;
ylim([min(summedLengthsPerFrame),max(summedLengthsPerFrame)]);

% Fitting
fit_idx = find(myTimes<FITTIME);
[fitMu A0] = DJK_ExponentialFit(myTimes(fit_idx)/60, summedLengthsPerFrame(fit_idx));

% Plotting
semilogy(myTimes,A0*2.^(fitMu*myTimes/60),'-k','LineWidth',2)
set(gca,'FontSize',15)

title([ 'Growth \newline' ...
        'Fitted \mu = ' num2str(fitMu) ' (dbl/hr)'])
xlabel('Time in minutes');
ylabel('Summed length all bacteria');

set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal')

saveas(h, [EXPORTFOLDER 'bacteriaSize_RAW-' identifier '.eps'])
saveas(h, [EXPORTFOLDER 'bacteriaSize_RAW-' identifier '.png'])

%% Fancy plot for suppl. mat. article
% Plotting setup
FONTSIZE = 24;

h = figure('pos', [100 100 700+100 600+100]); % 700x600
semilogy(myTimes, summedLengthsPerFrame,'o','LineWidth',4,'Color',[.6 .6 .6], 'MarkerSize',12);
hold on;
ylim([min(summedLengthsPerFrame)/2,max(summedLengthsPerFrame)*2]);

% Fitting
fit_idx = find(myTimes<FITTIME);
[fitMu A0] = DJK_ExponentialFit(myTimes(fit_idx)/60, summedLengthsPerFrame(fit_idx));

%title([ 'Colony growth'])
xlabel('Time (min)');
ylabel('Log summed lengths for colony (a.u.)');
set(gca,'YTickLabel',[]);
set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal')
set(gca,'FontSize',FONTSIZE)

saveas(h, [EXPORTFOLDER 'bacteriaSize-' identifier '.eps'])
saveas(h, [EXPORTFOLDER 'bacteriaSize-' identifier '.png'])

% Plotting of fitted line:
semilogy(myTimes,A0*2.^(fitMu*myTimes/60),'--k','LineWidth',5)

saveas(h, [EXPORTFOLDER 'bacteriaSize_withFit-' identifier '.eps'])
saveas(h, [EXPORTFOLDER 'bacteriaSize_withFit-' identifier '.png'])

disp('Done.');



%%

%{
%% Using Schnitzcells code

% copied from Excel file
p1 = DJK_initschnitz('Pos09-mini-01crop','2008-02-23','e.coli.amolf','rootDir','F:\X_Other_datasets\2008mondate\', 'cropLeftTop',[1 1], 'cropRightBottom',[800 800],'fluor1','none','fluor2','none','fluor3','none','setup','setup1','softwarePackage','metamorph','camera','?');
fitTime = DJK_analyzeMu(p1, schnitzcells, 'xlim', [0 2000], 'onScreen', 0,'fitTime',[1000 1200]);
%}









