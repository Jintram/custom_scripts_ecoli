
%%

if ~exist('PLOTOUTPUTDIRROOT','var')
    PLOTOUTPUTDIRROOT = 'D:\Local_Software\Martijn_extensions\Martijn_custom\filamentationrecovery\divisionmodel_plots\';
    % PLOTOUTPUTDIRROOT='D:\Local_Software\Martijn_extensions\Martijn_custom\filamentationrecovery\divisionmodel_plots\2016_10_05_conservativemodel\'
end

if ~exist('SAVEYESNO','var')
    SAVEYESNO=0; 
    % SAVEYESNO=1;
end

if ~exist('parameters','var')
    parameters.mu = 1.5; % dbls/hr
    parameters.replicationTime = 1/1.5; % hr/dbl
    parameters.replicationPhase = .1;
    parameters.divisionTime = 1; % for timer model to work, this needs to be equal to 1/mu
    parameters.divisionSize = 6;
    parameters.addedSize = 3;
    parameters.divisionType = 'nucleoidadder'; % simplesizer nucleoidsizer nucleoidadder
    parameters.ringType = 'nucleoid'; % nucleoid alwaysrandom alwaysmiddle alwaystypical dynamicalwaystypical
    parameters.divisionDelay = parameters.mu/10; % 0
    parameters.rechargeTime = parameters.mu/20; % 20 200
    parameters.divisionBlock=0;
    
    % noise parameters
    parameters.fluctuationIntensity = 0.04; % 0.04 seems good value
    parameters.relaxationTimeFluctuationsMinutes = 50; % 50 seems good value
    
    % Command to clear params:
    % clear parameters
end



INITIALLENGTH=3;
FIG3XLIM=[100 200];
SWITCHTIME = 150; % 150

SHOWNUCLEOIDS = 0;

HOWMANYSTARTINGSCHNITZES=100;

simulationtimes = 1:250; % simulation time

% extrapolated parametes
parameters.replicationTimeInMinutes = parameters.replicationTime*60;
parameters.lambda = parameters.mu*log(2); % /hr
parameters.lambdaperminute = parameters.lambda/60;
parameters.divisionTimePerMinute = parameters.divisionTime*60;
parameters.divisionDelayInMinutes = parameters.divisionDelay*60;
parameters.rechargeTimeInMinutes = parameters.rechargeTime*60;

% Initial conditions: create start of lineage
simulatedschnitzcells = {};
for schnitzIdx = 1:HOWMANYSTARTINGSCHNITZES
    simulatedschnitzcells(schnitzIdx).P = -1; % parent
    simulatedschnitzcells(schnitzIdx).D = 0; simulatedschnitzcells(schnitzIdx).E = 0; % daughters
    simulatedschnitzcells(schnitzIdx).cellLengths = [INITIALLENGTH];
    simulatedschnitzcells(schnitzIdx).growthEfficiencies = [1];
    simulatedschnitzcells(schnitzIdx).relativeRingSites = {[]};
    simulatedschnitzcells(schnitzIdx).times = [simulationtimes(1)];
    simulatedschnitzcells(schnitzIdx).nucleoidDuplicationEvents = [simulationtimes(1)- (1-parameters.replicationPhase).*parameters.replicationTimeInMinutes];
    simulatedschnitzcells(schnitzIdx).nrNucleoids = [1];
    simulatedschnitzcells(schnitzIdx).nucleoidDuplicationEventsBoolean=[1];
end
%% Establish output directory if desired

if SAVEYESNO==1

    currentDateTime = datestr(now);
    currentDateTime = strrep(currentDateTime, ':', '-');
    currentDateTime = strrep(currentDateTime, ' ', '_');
    plotoutputdir = [PLOTOUTPUTDIRROOT currentDateTime '\'];

    if ~exist(plotoutputdir)
        mkdir(plotoutputdir)
    end
end


%% simulation itself

% filamentation part of simulation
parameters.divisionBlock=1;
for time = simulationtimes(2:SWITCHTIME-1)
    [simulatedschnitzcells] = ODEdivisionlocation(simulatedschnitzcells, time, parameters);
    
    if mod(time,50)==0
        disp(['Simulated up to ' num2str(time) ' of ' num2str(simulationtimes(end)) '..']);
    end
end

% recovery part of simulation
parameters.divisionBlock=0;
for time = simulationtimes(SWITCHTIME:end)
    [simulatedschnitzcells] = ODEdivisionlocation(simulatedschnitzcells, time, parameters);
    
    if mod(time,50)==0
        disp(['Simulated up to ' num2str(time) ' of ' num2str(simulationtimes(end)) '..']);
    end
end


%%

h1=figure(1); clf; hold on;
plot([simulatedschnitzcells(:).times],[simulatedschnitzcells(:).cellLengths],'.');

plot([min(simulationtimes),max(simulationtimes)],[parameters.divisionSize,parameters.divisionSize],'-','Color',[.5 .5 .5])
plot([min(simulationtimes),max(simulationtimes)],[parameters.divisionSize,parameters.divisionSize]/2,'-','Color',[.5 .5 .5])

xlabel('Time [min]');
ylabel('Cell lengths [\mum]');

MW_makeplotlookbetter(20);

if SAVEYESNO
    saveas(h1,[plotoutputdir 'singlecellsize.fig']);
    saveas(h1,[plotoutputdir 'singlecellsize.tif']);
end

%% extra calculations

% calculate ring location & nucleoid location in absolute length values
for loopIdx = 1:numel(simulatedschnitzcells)   
    
    for idx = 1:numel(simulatedschnitzcells(loopIdx).times)
        
        % absolute ring sites
        simulatedschnitzcells(loopIdx).absoluteRingSite{idx}=...
            simulatedschnitzcells(loopIdx).cellLengths(idx).*simulatedschnitzcells(loopIdx).relativeRingSites{idx};
        
        % absolute nucleoid sites
        simulatedschnitzcells(loopIdx).absoluteNucleoidSite{idx}=...
            simulatedschnitzcells(loopIdx).cellLengths(idx).*(1:simulatedschnitzcells(loopIdx).nrNucleoids(idx))/(simulatedschnitzcells(loopIdx).nrNucleoids(idx)+1);
                
    end
    
end

%{
for time = simulationtimes
    
    indicesCurrentTime = find([simulatedschnitzcells(:).times] == time);
    
    allLengths = [simulatedschnitzcells(:).cellLengths];
    currentLengths = allLengths(indicesCurrentTime);
    totalLength = sum(currentLengths);
    cumsumLength = [0 cumsum(currentLengths)];
             
    %currentRingLocations=
    
%     plot(ones(numel(cumsumLength),1)*time, cumsumLength-totalLength/2,...
%         '-s','Color',[.5 .5 .5],...
%         'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[1,1,1],'LineWidth',4,'MarkerSize',4)
    
    plot(ones(numel(cumsumLength),1)*time, cumsumLength-totalLength/2,...
        'ok-','MarkerFaceColor','k','MarkerSize',5);
        %'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[1,1,1],'LineWidth',4)
    
end

xlabel('Time [min]');
ylabel('Cell lengths [\mum]');

MW_makeplotlookbetter(20);
%}

%%
% Create lookuptable
% lookuptable{n} gives, for nth simulation time, in lookuptable{n}(:,1) the
% schnitzes that are alive during that time, and lookuptable{n}(:,2), the
% corresponding frame in which they were alive. lookuptable{n}(:,3) is 0,
% but 1 if it was the last frame in which this schnitz was spotted.
%{
lookuptable={};
for t = 1:numel(simulationtimes)
    
    lookuptable{t} = [];
    hits=0;
    
    for i = 1:numel(simulatedschnitzcells)
    
        if any(simulatedschnitzcells(i).times==t)
            hits=hits+1;
            % schnitz#
            lookuptable{t}(hits,1) = [i];
            % which timefield corresponds
            lookuptable{t}(hits,2) = [find(simulatedschnitzcells(i).times==t)];
            % flag for if its the last frame in which schnitz lives
%             if (lookuptable{t}(hits,2))==numel(simulatedschnitzcells(i).times)
%                 lookuptable{t}(hits,3) = 1;
%             else
%                 lookuptable{t}(hits,3) = 0;
%             end
        end
    
    end
end
%}

%% build a lookuptable that corresponds with lineage structure
% lookuptable{n} gives, for nth simulation time, in lookuptable{n}(:,1) the
% schnitzes that are alive during that time, and lookuptable{n}(:,2), the
% corresponding number that corresponds to their progress towards their 
% life (i.e. schnitzcells(i).times==[time at n]). 
% 

lookuptable={}; lookuptable{1}=[]; hits=0;
for i = 1:numel(simulatedschnitzcells)
    if any(simulatedschnitzcells(i).times==simulationtimes(1))
        hits=hits+1;
        lookuptable{1}(hits,1) = i;
        lookuptable{1}(hits,2) = 1;
    end
end

tic;
% now expand lookuptable
for t = 2:numel(simulationtimes)
    
    lookuptable{t}=[]; 
    rprime = 0;
    for r = 1:numel(lookuptable{t-1}(:,1))
        
        % previous schnitz
        previousSchnitz = lookuptable{t-1}(r,1);
        
        % previous indices of time
        previousTimeIdx = find(simulatedschnitzcells(previousSchnitz).times==simulationtimes(t-1));
        
        % if previous time is last time
        if numel((simulatedschnitzcells(previousSchnitz).times))==previousTimeIdx
            % cell has divided
            nextSchnitz1 = simulatedschnitzcells(previousSchnitz).D;
            nextSchnitz2 = simulatedschnitzcells(previousSchnitz).E;
            % put in table
            rprime=rprime+1;
            lookuptable{t}(rprime,1) = nextSchnitz1;
            lookuptable{t}(rprime,2) = 1;
            rprime=rprime+1;
            lookuptable{t}(rprime,1) = nextSchnitz2;
            lookuptable{t}(rprime,2) = 1;
        else
            % cell has not divided
            nextSchnitz = previousSchnitz;
            rprime=rprime+1;
            lookuptable{t}(rprime,1) = nextSchnitz;
            lookuptable{t}(rprime,2) = find((simulatedschnitzcells(nextSchnitz).times==simulationtimes(t)));
        end
            
        if toc>10
            error('timeout error');
        end
        
    end
    
end

%%

h3=figure(3); clf; hold on;
% dummy markers for legend (appear at t<0)
l1=plot(-1,0,'sk-','MarkerFaceColor','k','MarkerSize',5);
l2=plot(-1,0,'or','MarkerSize',4,'MarkerFaceColor','r');
l3=plot(-1,0,'ob','MarkerSize',4);

tic
for t=1:numel(lookuptable)
    
    currentSchnitzes = lookuptable{t}(:,1)';
    currentCellPhaseIdx = lookuptable{t}(:,2)';
    
    totalLength=0;
    for loopIdx = 1:numel(currentSchnitzes)
        totalLength = totalLength+simulatedschnitzcells(currentSchnitzes(loopIdx)).cellLengths(currentCellPhaseIdx(loopIdx));
    end
    
    previousLengthsSummed=0;
    plottingDividedLocations = [-totalLength/2];
    plottingRingLocations = []; plottingNucleoidLocations= [];
    for loopIdx = 1:numel(currentSchnitzes)
    
        plottingDividedLocations = ...
            [plottingDividedLocations previousLengthsSummed+simulatedschnitzcells(currentSchnitzes(loopIdx)).cellLengths(currentCellPhaseIdx(loopIdx))-totalLength/2];
        
        for ringIdx = 1:numel(simulatedschnitzcells(currentSchnitzes(loopIdx)).absoluteRingSite{currentCellPhaseIdx(loopIdx)})
            plottingRingLocations = ...
                [plottingRingLocations previousLengthsSummed+simulatedschnitzcells(currentSchnitzes(loopIdx)).absoluteRingSite{currentCellPhaseIdx(loopIdx)}(ringIdx)-totalLength/2];
        end

        for nucleoidIdx = 1:numel(simulatedschnitzcells(currentSchnitzes(loopIdx)).absoluteNucleoidSite{currentCellPhaseIdx(loopIdx)})
            plottingNucleoidLocations = ...
                [plottingNucleoidLocations previousLengthsSummed+simulatedschnitzcells(currentSchnitzes(loopIdx)).absoluteNucleoidSite{currentCellPhaseIdx(loopIdx)}(nucleoidIdx)-totalLength/2];
        end
        
        %plot(simulatedschnitzcells(currentSchnitzes(loopIdx)).times(currentFrames(loopIdx)),...
        %     previousLengthsSummed+simulatedschnitzcells(currentSchnitzes(loopIdx)).cellLengths(currentFrames(loopIdx))-totalLength/2,...
        %    'ok-','MarkerFaceColor','k','MarkerSize',5);
        
        previousLengthsSummed=previousLengthsSummed+simulatedschnitzcells(currentSchnitzes(loopIdx)).cellLengths(currentCellPhaseIdx(loopIdx));
       
        if toc>10
            error('timeout error');
        end
        
    end
    
    plot(ones(1,numel(plottingDividedLocations))*simulationtimes(t),...
            plottingDividedLocations,...
           'sk-','MarkerFaceColor','k','MarkerSize',5);
    plot(ones(1,numel(plottingRingLocations))*simulationtimes(t),...
            plottingRingLocations,...
           'or','MarkerSize',4,'MarkerFaceColor','r');
    if SHOWNUCLEOIDS
        plot(ones(1,numel(plottingNucleoidLocations))*simulationtimes(t),...
                plottingNucleoidLocations,...
               'ob','MarkerSize',4);
    end
end

currentYlim=ylim;
plot([SWITCHTIME,SWITCHTIME], currentYlim,'--k');

xlim(FIG3XLIM)

xlabel('Time');
ylabel('Cells [um]');

if SHOWNUCLEOIDS
    legend([l1,l2,l3],{'Cells','Division ring','Nucleoid'},'Location','SouthWest'); % SouthWest Best
else
    legend([l1,l2],{'Cells','Division ring'},'Location','SouthWest'); % SouthWest Best
end
warning('This does not always work because lines are not always plotted during above loop..');

MW_makeplotlookbetter(20);
%%
if SAVEYESNO
    saveas(h3,[plotoutputdir 'cellview.fig']);
    saveas(h3,[plotoutputdir 'cellview.tif']);
end

%% histogram of sizes

h4=figure(4); clf; hold on;

allLengths = [simulatedschnitzcells(:).cellLengths];
[counts, locations]=hist(allLengths,100);
meanLength = mean(allLengths);
medianLength = median(allLengths);
plot(locations, counts);
plot([medianLength,medianLength],[0,max(counts)],'-r')
title(['median = ' num2str(medianLength)]);

if SAVEYESNO
    saveas(h4,[plotoutputdir 'cellsizes.fig']);
    saveas(h4,[plotoutputdir 'cellsizes.tif']);
end

%% Save parameter struct
if SAVEYESNO
    save([plotoutputdir 'parameters.mat'],'parameters');
end

% output parameters to text file
stringwithdata = '';
allFieldNames = {fieldnames(parameters)};
for i = 1:numel(allFieldNames{1})
    fieldName = allFieldNames{1}{i}
    fieldValue = parameters.(fieldName);
    stringwithdata = [stringwithdata 10 fieldName '=' num2str(fieldValue)];
end
stringwithdata
%%
if SAVEYESNO
    fid = fopen([plotoutputdir 'parameters.txt'],'w');
    fprintf(fid, stringwithdata);
    fclose(fid);
end


%{

schnitzcells=simulatedschnitzcells;
save('D:\Local_Software\Martijn_extensions\Martijn_custom\filamentationrecovery\simulatedSchnitzcells\schnitz1_longrecharge.mat','schnitzcells','parameters');

%}







