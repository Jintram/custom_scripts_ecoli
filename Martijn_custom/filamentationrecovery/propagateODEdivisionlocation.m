
%%
parameters.mu = 1;
parameters.divisionTime = 1; % for timer model to work, this needs to be equal to mu
parameters.divisionSize = 6;
parameters.addedSize = 6;
parameters.divisionType = 'adder';

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
h=figure(3); clf; hold on;
% Create lookuptable
% lookuptable{n} gives, for nth simulation time, in lookuptable{n}(:,1) the
% schnitzes that are alive during that time, and lookuptable{n}(:,2), the
% corresponding frame in which they were alive.
lookuptable={};
lineageHierarchy = {};
for t = simulationtimes    
    
    lookuptable{t} = [];
    hits=0;
    
    for i = 1:numel(simulatedschnitzcells)
    
        if any(simulatedschnitzcells(i).times==t)
            hits=hits+1;
            lookuptable{t}(hits,1) = [i];
            lookuptable{t}(hits,2) = [find(simulatedschnitzcells(i).times==t)];
        end
    
    end
end

% now sort the lookuptable by lineage hierarchy
% go over times again
for t = simulationtimes(2:end)
    
    % go over schnitzes
    for tableRow = 1:numel(lookuptable{t}(:,1))
        
        
        
        % if it is newborn, move it to the location in the array it's
        % parent was previously
        if lookuptable{t}(tableRow,2)==1 % if time was 1st element
            
            parentLocation = find(lookuptable{t-1}(tableRow,1)==simulatedschnitzcells(tableRow).P);
            originalLocation = tableRow;
            
            % remember old values
            oldValues = lookuptable{t}(tableRow,:);
            % remove at original location
            lookuptable{t}=lookuptable{t}([1:tableRow-1,tableRow+1:end],:);
            % insert at new location
            lookuptable{t}=[lookuptable{t}([1:tableRow],:);oldValues;lookuptable{t}([tableRow:end],:)];            
            
            
        end
    end
    
end

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
