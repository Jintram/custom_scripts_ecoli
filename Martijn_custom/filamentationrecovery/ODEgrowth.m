



function [cellLengths,divisionEvents,divisionSizes,growthEfficiencies] = ODEgrowth(cellLengths,time,divisionEvents,divisionSizes,growthEfficiencies,parameters)

%% parameters

if ~isfield(parameters,'mu')
    parameters.mu = 1;    
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
    
divide=0;
    
%% cell division

% dt = 1 minute, and is implicit

% sizer
if strcmp(parameters.divisionType,'sizer')
    if cellLengths(end)>parameters.divisionSize
        divide = 1;   
    end
end

% adder
if strcmp(parameters.divisionType,'adder')
    if cellLengths(end)>parameters.divisionSize
        lastDivisionSize = divisionSizes(end);
        if (cellLengths(end)-lastDivisionSize)>parameters.addedSize
            divide=1;
        end
    end
end

% timer
if strcmp(parameters.divisionType,'timer')
    if (time-divisionEvents(end))>parameters.divisionTimePerMinute
        divide=1;
    end
end

% execute division if determined so
if divide
    cellLengths(end)=cellLengths(end)/2;
    
    divisionEvents(end+1) = time;
    divisionSizes(end+1) = cellLengths(end);    
end

%% ODEs
deltaGrowthEfficiency = -(growthEfficiencies(end)-1)*1/parameters.relaxationTimeFluctuationsMinutes + ... % damping; 1 is the equilibrium value
                        ((rand()-.5)*2)*parameters.fluctuationIntensity; % noise component

deltacellLength = growthEfficiencies(end)*parameters.lambdaperminute*cellLengths(end);                    
                    
% propagate
growthEfficiencies(end+1) = growthEfficiencies(end)+deltaGrowthEfficiency;
cellLengths(end+1) = cellLengths(end) + deltacellLength;




