
%%
parameters.mu = 1;
parameters.divisionTime = 1; % for timer model to work, this needs to be equal to mu
parameters.divisionSize = 6;
parameters.addedSize = 6;
parameters.divisionType = 'sizer';

% noise parameters
parameters.fluctuationIntensity = 0.04; % 0.04 seems good value
parameters.relaxationTimeFluctuationsMinutes = 50; % 50 seems good value

INITIALLENGTH=25;

simulationtimes = 1:100; % simulation time

% Create start of lineage
simulatedschnitzcells = {};
simulatedschnitzcells(1).P = -1; % parent
simulatedschnitzcells(1).D = 0; simulatedschnitzcells(1).E = 0; % daughters
simulatedschnitzcells(1).cellLengths = [INITIALLENGTH];
simulatedschnitzcells(1).growthEfficiencies = [1];
simulatedschnitzcells(1).relativeRingSites = {[.5]};
simulatedschnitzcells(1).times = [simulationtimes(1)];

%%
for time = simulationtimes(2:end)
    [simulatedschnitzcells] = ODEdivisionlocation(simulatedschnitzcells, time, parameters);
    
    if mod(time,50)==0
        disp(['Simulated up to ' num2str(time) ' of ' num2str(simulationtimes(end)) '..']);
    end
end

%%

h=figure(1); clf; hold on;
plot([simulatedschnitzcells(:).times],[simulatedschnitzcells(:).cellLengths],'.');

plot([min(simulationtimes),max(simulationtimes)],[parameters.divisionSize,parameters.divisionSize],'-','Color',[.5 .5 .5])
plot([min(simulationtimes),max(simulationtimes)],[parameters.divisionSize,parameters.divisionSize]/2,'-','Color',[.5 .5 .5])


%% kymograph
h=figure(2); clf; hold on;

% calculate ring location in absolute length values
for loopIdx = 1:numel(simulatedschnitzcells)   
    
    for idx = 1:numel(simulatedschnitzcells(loopIdx).times)
        simulatedschnitzcells(loopIdx).absoluteRingSite{idx}=...
            simulatedschnitzcells(loopIdx).cellLengths(idx).*simulatedschnitzcells(loopIdx).relativeRingSites{idx};
    end
    
end

for time = simulationtimes(1:100)
    
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
% corresponding frame in which they were alive. lookuptable{n}(:,3) is 0,
% but 1 if it was the last frame in which this schnitz was spotted.
% identify the first schnitzcells
lookuptable={}; lookuptable{1}=[]; hits=0;
for i = 1:numel(simulatedschnitzcells)
    if any(simulatedschnitzcells(i).times==simulationtimes(1))
        hits=hits+1;
        lookuptable{1}(hits,1) = i;
        lookuptable{1}(hits,2) = 1;
    end
end

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
            
    end
    
end

%%
h=figure(3); clf; hold on;
for t=1:numel(lookuptable)
    
    currentSchnitzes = lookuptable{t}(:,1)';
    currentFrames    = lookuptable{t}(:,2)';
    
    totalLength=0;
    for loopIdx = 1:numel(currentSchnitzes)
        totalLength = totalLength+simulatedschnitzcells(currentSchnitzes(loopIdx)).cellLengths(currentFrames(loopIdx));
    end
    
    previousLengthsSummed=0;
    plottingDividedLocations = [-totalLength/2];
    plottingRingLocations = [];
    for loopIdx = 1:numel(currentSchnitzes)
    
        plottingDividedLocations = ...
            [plottingDividedLocations previousLengthsSummed+simulatedschnitzcells(currentSchnitzes(loopIdx)).cellLengths(currentFrames(loopIdx))-totalLength/2];
        
        for ringIdx = 1:numel(simulatedschnitzcells(currentSchnitzes(loopIdx)).absoluteRingSite{currentFrames(loopIdx)})
            plottingRingLocations = ...
                [plottingRingLocations previousLengthsSummed+simulatedschnitzcells(currentSchnitzes(loopIdx)).absoluteRingSite{currentFrames(loopIdx)}(ringIdx)-totalLength/2];
        end
        
        %plot(simulatedschnitzcells(currentSchnitzes(loopIdx)).times(currentFrames(loopIdx)),...
        %     previousLengthsSummed+simulatedschnitzcells(currentSchnitzes(loopIdx)).cellLengths(currentFrames(loopIdx))-totalLength/2,...
        %    'ok-','MarkerFaceColor','k','MarkerSize',5);
        
        previousLengthsSummed=previousLengthsSummed+simulatedschnitzcells(currentSchnitzes(loopIdx)).cellLengths(currentFrames(loopIdx));
        
    end
    
    plot(ones(1,numel(plottingDividedLocations))*simulationtimes(t),...
            plottingDividedLocations,...
           'sk-','MarkerFaceColor','k','MarkerSize',5);
    plot(ones(1,numel(plottingRingLocations))*simulationtimes(t),...
            plottingRingLocations,...
           'or','MarkerSize',4);
end


%{
figure(1); clf; hold on;
plot(times,cellLengths,'LineWidth',3)

plot(divisionEvents,zeros(1,numel(divisionEvents)),'^r','MarkerFaceColor','r','MarkerSize',15)

xlabel('time (min)');
ylabel('cellength [\mum]');
MW_makeplotlookbetter(20);

figure(2); clf; hold on;
plot(times,growthEfficiencies,'LineWidth',3)

plot(divisionEvents,zeros(1,numel(divisionEvents)),'^r','MarkerFaceColor','r','MarkerSize',15)

xlabel('time (min)');
ylabel('growth efficiency');
MW_makeplotlookbetter(20);
%}
