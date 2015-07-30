function X = SSA_direct_function(stoichiometry,prop,speciesCounts0,probeTimeArray,Nspecies)
%SSA_direct_function Calculates gillespie simulation.
%   SSA_direct_function gives the trajectory X of reactions 
%
% Written by Johannes Keegstra, adpated by Martijn Wehrens
% 2015/02

% Prepare output vector
X = NaN(Nspecies,length(probeTimeArray));

% Prepare initial state
time = 0;
speciesCounts = speciesCounts0;
idx3 = 1;
X(:,idx3) = speciesCounts0;

% For message
reverseStr='';
% Performance
maxprobeTimeArray = max(probeTimeArray);
totalSteps = numel(probeTimeArray); % TODO
% Track performance
totalTimePassed = 0; tic;
% Simulation loop
while time < maxprobeTimeArray
    
    % Calculate propensity matrix
    propensityMatrix = prop(speciesCounts);
    
    % Get time to first rxn by assuming k=sum_i(k_i), exponential fn
    % Determine tau stochastically
    taumin = -log(rand(1,1))/sum(propensityMatrix); 
    
    % Determine which reaction fires stochastically
    randomnr = rand(1,1)*sum(propensityMatrix); 
    rxnToFire=1;
    while sum(propensityMatrix(1:rxnToFire)) < randomnr
        rxnToFire=rxnToFire+1;
    end
    
    % If we landed in the next step, record and proceed
    if time > probeTimeArray(idx3);
        % track performance
        timepassed = toc;
        
        % actual algorithm
        idx3=idx3+1;
        X(:,idx3)=speciesCounts;

       % print progress to user
       totalTimePassed = totalTimePassed + timepassed;
       ETA = (totalTimePassed/idx3) * (totalSteps-idx3);
       msg = sprintf('Processed time %.2f/%.2f s or %d/%d steps in %.2f minutes (ETA %.2f mins)', time, maxprobeTimeArray, idx3, totalSteps, totalTimePassed/60, ETA/60);
       fprintf([reverseStr, msg]);
       reverseStr = repmat(sprintf('\b'), 1, length(msg));
       
       % track performance
       tic;
    end
        
    % Progress time
    time=time+taumin;
    
    % If end of simulation not reached, update speciescounts according 
    % sto
    if time>(max(probeTimeArray));
          break
    else
           speciesCounts=speciesCounts+stoichiometry(:,rxnToFire); %update reactions.
    end
    
  
end



end

