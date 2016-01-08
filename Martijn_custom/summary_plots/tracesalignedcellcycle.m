% I was wondering whether there is a systematic effect during the cell
% cycle visible in growth. Note that Noreen Walker has also looked into
% this.
% Here I just make some plots of cellular traces with t0 = division.

%% Load some schnitzcells dataset
if ~exist('MYDATASET','var')
    MYDATASET='F:\A_Tans1_step1_incoming_not_backed_up\2015-10-20\pos4crop\data\pos4crop-Schnitz.mat';
    disp('Loaded default dataset');
end
load(MYDATASET);

if ~exist('schnitzcells','var')    
    schnitzcells=s_rm;
end

if ~exist('speedfield')
    speedfield = 'mu5_fitNew';
end

%% Get all cell traces into one cell (this is actually not necessary as we
% could have also used "schnitzcells.mu5_fitNew" all the time).
traces = {schnitzcells.(speedfield)};
maxTraceLength = max(cellfun('length',traces)); % max trace length in frames

% make figure
figure(1); clf; hold on;

%% group data per frame
musPerFrame = cell(1,maxTraceLength);  
extendedtraces = [];
% extendedtraces = nan(numel(traces),maxTraceLength); % not convenient
% since we need to skip empty ones
for i = 1:numel(traces)
    
    if ~isempty(traces{i})
        extendedtraces = [extendedtraces; imresize(traces{i},[1,maxTraceLength],'nearest')];
    end
    
    for j = 1:numel(traces{i})
        musPerFrame{j} = [musPerFrame{j} traces{i}(j)];
    end
    plot(traces{i},'-')
end

%% violin plots
figure(2); clf; hold on;
%violin(musPerFrame);
violin(extendedtraces);

%% calculate summary quantity
stdovermeanPerFrame = [];
stdovermeanPerFrameForExtended = [];
for i=1:maxTraceLength
    stdovermeanPerFrame(i) = std(musPerFrame{i})/mean(musPerFrame{i});
    indexNoNaN=~isnan(extendedtraces(1,:));
    stdovermeanPerFrameForExtended(i) = std(extendedtraces(indexNoNaN,i))/mean(extendedtraces(indexNoNaN,i));
end

% plot of summary parameter
figure(3); clf; hold on;
%h=bar(stdovermeanPerFrame)
h=bar(stdovermeanPerFrameForExtended);
set(h, 'FaceColor', [.7 .7 .7]);

ylabel('Coefficient of variation (CV) (std/mean)');
xlabel('Time in cell cycle [frames]');
MW_makeplotlookbetter(15);

xlim([0, maxTraceLength+1]);
ylim([0,max(stdovermeanPerFrameForExtended)*1.1]);



