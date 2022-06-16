function [simulatedschnitzcells] = ODEdivisionlocation(simulatedschnitzcells, time, parameters)

%% Description
%
% This function simulates a lineage, focusing on division locations.
%
% Input and output is "simulatedschnitzcells", which has a similar
% structure to the schnitzcells resulting from experimental analysis.
%
%
% - MW 2016/06


%% parameters

% overflow prevention errors
MAXNUCLEOIDSINCELL = 1000;
MAXSCHNITZESTOMAKE = 10000;

    
    
%% Go over active schnitzcells

% first identify them
activeSchnitzcells = find( ([simulatedschnitzcells(:).D] == 0) & ([simulatedschnitzcells(:).E] == 0) );

for schnitzIdx = activeSchnitzcells
      
    %% ODEs ---------------------------------------------------------------
    deltaGrowthEfficiency = -(simulatedschnitzcells(schnitzIdx).growthEfficiencies(end)-1)*1/parameters.relaxationTimeFluctuationsMinutes + ... % damping; 1 is the equilibrium value
                            ((rand()-.5)*2)*parameters.fluctuationIntensity; % noise component

    deltacellLength = simulatedschnitzcells(schnitzIdx).growthEfficiencies(end)*parameters.lambdaperminute*simulatedschnitzcells(schnitzIdx).cellLengths(end);                    

    
    % prepare normal propagation
    newtime = ...
        time;        
    newGrowthEfficiency = ...
        simulatedschnitzcells(schnitzIdx).growthEfficiencies(end)+deltaGrowthEfficiency;
    newCellLength        = ...
        simulatedschnitzcells(schnitzIdx).cellLengths(end) + deltacellLength;
    
    % generally, nr. nucleoids and rings are identical to last frame
    newNrNucleoids = ...
        simulatedschnitzcells(schnitzIdx).nrNucleoids(end);
    newRingSites = ...
        simulatedschnitzcells(schnitzIdx).relativeRingSites{end};
    
    % determine last nucleoid duplication event
    nucleoidDuplicationEvents=[]; checkoutSchnitz = schnitzIdx;
    while isempty(nucleoidDuplicationEvents)
        nucleoidDuplicationEvents = simulatedschnitzcells(checkoutSchnitz).nucleoidDuplicationEvents;
        checkoutSchnitz = simulatedschnitzcells(checkoutSchnitz).P;
    end
    lastDuplicationEvent=nucleoidDuplicationEvents(end);
    
    % determine last division event
    lastDivisionEvent = simulatedschnitzcells(schnitzIdx).times(1);
    
    
    %% determine whether nucleoid duplication & cell division should happen
    % dt = 1 minute, and is implicit
    
    % Do not duplicate nucleoids, unless..
    nucleoidDuplication=0;
    % Do not divide, unless..
    divide=0;
           
    % Determine according to divison type of the cell
    if strcmp(parameters.divisionType,'nucleoidsizer') 
    
        % nucleoid sizer
        % ===
    
        % nucleoid duplication
        if newCellLength>(parameters.divisionSize)  && ... % sizer condition
            (newtime-(lastDuplicationEvent)) > (parameters.replicationTimeInMinutes) % timer condition
        
            nucleoidDuplication=1;
            lastDuplicationEvent=newtime;
        end
        if parameters.divisionDelayInMinutes< (newtime-lastDuplicationEvent) ... % delay after duplication
                && parameters.rechargeTimeInMinutes< (newtime-lastDivisionEvent) ... delay after previous division
                && ~parameters.divisionBlock
        % cell division
        %if newCellLength>parameters.divisionSize
            divide = 1;   
        end    
    
    elseif strcmp(parameters.divisionType,'nucleoidadder')
    
        % nucleoid adder
        % ===
        
        lastDivisionSize = simulatedschnitzcells(schnitzIdx).cellLengths(1);
        
        % nucleoid duplication
        %if newCellLength>parameters.divisionSize
        if (newCellLength-lastDivisionSize)>(parameters.addedSize) && ... % adder condition
            (newtime-(lastDuplicationEvent)) > (parameters.replicationTimeInMinutes) % timer condition
        
            nucleoidDuplication=1;
            lastDuplicationEvent=newtime;
        end
        %end
        
        % cell division
        %if newCellLength>parameters.divisionSize
        if parameters.divisionDelayInMinutes< (newtime-lastDuplicationEvent)... % delay after duplication
                && parameters.rechargeTimeInMinutes< (newtime-lastDivisionEvent) ... delay after previous division
                && ~parameters.divisionBlock
            if (newCellLength-lastDivisionSize)>parameters.addedSize
                divide=1;
            end
        end
        
    elseif strcmp(parameters.divisionType,'nucleoidtimer')
    
        % nucleoid timer
        % ===
    
        % nucleoid duplication
        if (newtime-(lastDuplicationEvent)) > (parameters.replicationTimeInMinutes) % timer condition
            
            nucleoidDuplication=1;
            lastDuplicationEvent=newtime;
        end
        % division
        %if (simulatedschnitzcells(schnitzIdx).times(1)-newtime)>parameters.divisionTimePerMinute
        if parameters.divisionDelayInMinutes< (newtime-lastDuplicationEvent) ... % delay after duplication
                && parameters.rechargeTimeInMinutes< (newtime-lastDivisionEvent) ... delay after previous division
                && ~parameters.divisionBlock
            divide=1;
        end
        
    elseif strcmp(parameters.divisionType,'simplesizer')        
        
        if newCellLength>(parameters.divisionSize) ... % simple sizer condition
                && parameters.rechargeTimeInMinutes< (newtime-lastDivisionEvent) ... % delay after previous division
                && ~parameters.divisionBlock % cell division inhibition
            divide=1;
        end
        
    end    
    
    

    %% execute duplication event if determined so
        % nucleoid replication that results in more rings (assuming all
    % nucleoids divide in phase and are thus automatically aligned with 
    % rings).        
    if nucleoidDuplication
               
        % duplicate nucleoid number
        newNrNucleoids = newNrNucleoids*2;        
        
        % update duplication event list
        simulatedschnitzcells(schnitzIdx).nucleoidDuplicationEvents(end+1) = newtime;        
    end
    
    %% Update rings accordingly -------------------------------------------
    
    if strcmp(parameters.ringType,'nucleoid')
        
        % rings follow nucleoid number
        if nucleoidDuplication            
            nrRings = newNrNucleoids-1;        
            newRingSites = (1:nrRings)/(nrRings+1);
        end
    
    elseif strcmp(parameters.ringType,'alwaysmiddle')
        
        if divide
            newRingSites = .5;
        end
        
    elseif strcmp(parameters.ringType,'alwaystypical')
        
        % rings are such that resulting length closest to (multiple) typical
        nrRings = round(newCellLength/(parameters.divisionSize/2))-1;
        newRingSites = (1:nrRings)/(nrRings+1);
    
    elseif strcmp(parameters.ringType,'dynamicalwaystypical')
        
        % rings are such that resulting length closest to (multiple) typical
        nrRings = round(newCellLength/(parameters.divisionSize/2))-1;
        newRingSites = (1:nrRings)/(nrRings+1);
        newRingSites = newRingSites(find(round(rand(numel(newRingSites),1)))); % select half of them to simulate dynamic behavior
        
    elseif strcmp(parameters.ringType,'alwaysrandom')
        
        marginAtSides = (parameters.divisionSize/4 * .99) / newCellLength;
        widthOfPotentialRingSites = 1-2*marginAtSides;
        newRingSites = rand()*widthOfPotentialRingSites+marginAtSides;
        
    end
    
    %% execute division if determined so ----------------------------------
    % note that life for a cell starts here, so the time is recorded in the
    % new cells, but not in the old ones.
    divided=0;
    if divide
        
        % choose (at random) at which ring to divide
        nrRings = numel(newRingSites);
        
        % divide if there's a ring
        if nrRings > 0
            selectedRing = ceil(rand*nrRings);
            ratioWhereToDivide = newRingSites(selectedRing);

            % create daughter indices
            newIdx1 = numel(simulatedschnitzcells)+1;
            newIdx2 = numel(simulatedschnitzcells)+2;

            % link mother and daughters
            simulatedschnitzcells(schnitzIdx).D=newIdx1;
            simulatedschnitzcells(schnitzIdx).E=newIdx2;
            simulatedschnitzcells(newIdx1).P=schnitzIdx;
            simulatedschnitzcells(newIdx2).P=schnitzIdx;

            % New cells don't have children
            simulatedschnitzcells(newIdx1).D=0;
            simulatedschnitzcells(newIdx1).E=0;
            simulatedschnitzcells(newIdx2).D=0;
            simulatedschnitzcells(newIdx2).E=0;

            % calculate daughter cell lengths
            lengthBeforeDivision = newCellLength;
            lengthAfterDivision1 = lengthBeforeDivision*ratioWhereToDivide;
            lengthAfterDivision2 = lengthBeforeDivision-lengthAfterDivision1; %lengthBeforeDivision*(1-ratioWhereToDivide);
            simulatedschnitzcells(newIdx1).cellLengths=lengthAfterDivision1;
            simulatedschnitzcells(newIdx2).cellLengths=lengthAfterDivision2;
            % daughter growthEfficiencies
            simulatedschnitzcells(newIdx1).growthEfficiencies=[newGrowthEfficiency];
            simulatedschnitzcells(newIdx2).growthEfficiencies=[newGrowthEfficiency];
            % daughter times
            simulatedschnitzcells(newIdx1).times=time;
            simulatedschnitzcells(newIdx2).times=time;     
            % daughter ring reorganization times
            simulatedschnitzcells(newIdx1).ringReorganizationTimes = [];
            simulatedschnitzcells(newIdx1).ringReorganizationTimes = [];
            % daughter nr of nucleoids (assume distributed by length)
            nrNucleoidsForCell1 = round(newNrNucleoids*lengthAfterDivision1/lengthBeforeDivision);
            simulatedschnitzcells(newIdx1).nrNucleoids = nrNucleoidsForCell1;
            simulatedschnitzcells(newIdx2).nrNucleoids = newNrNucleoids-nrNucleoidsForCell1;
            % nucleoid division boolean event list
            simulatedschnitzcells(newIdx1).nucleoidDuplicationEventsBoolean(end+1) = nucleoidDuplication;
            simulatedschnitzcells(newIdx2).nucleoidDuplicationEventsBoolean(end+1) = nucleoidDuplication;
            
            
            % daughter ring sites        
            % (compensate for )
            simulatedschnitzcells(newIdx1).relativeRingSites = ...
                { newRingSites(1:selectedRing-1)   .* lengthBeforeDivision / lengthAfterDivision1};
            simulatedschnitzcells(newIdx2).relativeRingSites = ...
                { ((newRingSites(selectedRing+1:end) .* lengthBeforeDivision) - lengthAfterDivision1) / lengthAfterDivision2};

            divided=1;
        else
            disp('No ring present to divide..');
        end
        
    end
        
    if ~divided 
        
        % propagate normally
        simulatedschnitzcells(schnitzIdx).times(end+1) = newtime;
        simulatedschnitzcells(schnitzIdx).growthEfficiencies(end+1) = ...
            newGrowthEfficiency;
        simulatedschnitzcells(schnitzIdx).cellLengths(end+1)        = ...
            newCellLength;   
        
        simulatedschnitzcells(schnitzIdx).relativeRingSites{end+1} = ...
            newRingSites;
        
        simulatedschnitzcells(schnitzIdx).nrNucleoids(end+1) = ...
            newNrNucleoids;

        simulatedschnitzcells(schnitzIdx).nucleoidDuplicationEventsBoolean(end+1) = ...
            nucleoidDuplication;
        
    end
    
    if newNrNucleoids > MAXNUCLEOIDSINCELL
        error('Overflow protection error (nucleoids).');
    end
    
end

if numel(simulatedschnitzcells)>MAXSCHNITZESTOMAKE
    error('Overflow protection error (schnitznr).');
end


end