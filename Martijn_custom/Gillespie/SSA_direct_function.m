function X = SSA_direct_function(stoichiometry,prop,speciesCounts0,probeTimeArray)
%SSA_direct_function Calculates gillespie simulation.
%   SSA_direct_function gives the trajectory X of reactions 
%
% Written by Johannes Keegstra, adpated by Martijn Wehrens
% 2015/02

% Prepare output vector
X = NaN(2,length(probeTimeArray));

% Prepare initial state
time = 0;
speciesCounts = speciesCounts0;
idx3 = 1;
X(:,idx3) = speciesCounts0;

% Simulation loop
while time < (max(probeTimeArray))
    
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
        idx3=idx3+1;
        X(:,idx3)=speciesCounts;
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

