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

if ~isfield(parameters,'mu')
    parameters.mu = 1;    
end
if ~isfield(parameters,'divisionTimePerMinute')
    parameters.divisionTime = 1;    
end
if ~isfield(parameters,'divisionSize')
    parameters.divisionSize = 6;    
end
if ~isfield(parameters,'divisionType')
    parameters.divisionType = 'sizer';    
end
if ~isfield(parameters,'fluctuationIntensity')
    parameters.fluctuationIntensity = 0;    
end
if ~isfield(parameters,'relaxationTimeFluctuationsMinutes')
    parameters.relaxationTimeFluctuationsMinutes = 1;    
end

parameters.lambda=parameters.mu*log(2);
parameters.lambdaperminute = parameters.lambda/60;
parameters.divisionTimePerMinute = parameters.divisionTime*60;
    
    
%% Go over active schnitzcells

% first identify them
activeSchnitzcells = find( ([simulatedschnitzcells(:).D] == 0) & ([simulatedschnitzcells(:).E] == 0) );

for schnitzIdx = activeSchnitzcells
      
    %% ODEs
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
    
    % determine potential division sites
    newRingSites = ...
        simulatedschnitzcells(schnitzIdx).relativeRingSites{end};
    
    %% cell division
    % dt = 1 minute, and is implicit
    
    % Do not divide, unless..
    divide=0;
       
    % sizer
    if strcmp(parameters.divisionType,'sizer')
        if newCellLength>parameters.divisionSize
            divide = 1;   
        end
    end

    % adder
    if strcmp(parameters.divisionType,'adder')
        if newCellLength>parameters.divisionSize
            lastDivisionSize = simulatedschnitzcells(schnitzIdx).cellLengths(1);
            if (newCellLength-lastDivisionSize)>parameters.addedSize
                divide=1;
            end
        end
    end

    % timer
    if strcmp(parameters.divisionType,'timer')
        if (time-newtime)>parameters.divisionTimePerMinute
            divide=1;
        end
    end

    % execute division if determined so
    if divide
        
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
        simulatedschnitzcells(newIdx1).cellLengths=[newCellLength*newRingSites];
        simulatedschnitzcells(newIdx2).cellLengths=[newCellLength*(1-newRingSites)];
        % daughter growthEfficiencies
        simulatedschnitzcells(newIdx1).growthEfficiencies=[newGrowthEfficiency];
        simulatedschnitzcells(newIdx2).growthEfficiencies=[newGrowthEfficiency];
        % daughter times
        simulatedschnitzcells(newIdx1).times=time;
        simulatedschnitzcells(newIdx2).times=time;            
    
        % daughter ring sites        
        simulatedschnitzcells(newIdx1).relativeRingSites = ...
            {newRingSites};
        simulatedschnitzcells(newIdx2).relativeRingSites = ...
            {newRingSites};
        
    else       
        
        % propagate normally
        simulatedschnitzcells(schnitzIdx).times(end+1) = newtime;
        simulatedschnitzcells(schnitzIdx).growthEfficiencies(end+1) = ...
            newGrowthEfficiency;
        simulatedschnitzcells(schnitzIdx).cellLengths(end+1)        = ...
            newCellLength;   
        
        simulatedschnitzcells(schnitzIdx).relativeRingSites{end+1} = ...
            newRingSites;
        
    end
    
end


end