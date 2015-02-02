
% Written by Johannes Keegstra, adapted by Martijn Wehrens
% 2015/02

clear all

% Specify initial and final times
time    = 0;
tstop   = 10; 

% Initial concentrations of species
speciesCounts0      = [10; 0]; 

% Stoichiometry matrix (e.g. 2A -> B)
stoichiometry = [1 -1 0 0; ...
                 0 0 1 -1]; 

% Define times at which simulation is halted to probe concentrations 
probeNTimes = 10;
probeTimeArray = linspace(0,tstop,probeNTimes);  % times at which to plot.

% Calculate propensities
kmrnaProd = 10; kmrnaDeg = 1; kProteinProd = 10; kProteinDeg = 1;
prop = @(x)([kmrnaProd;             % mRNA prod. rate
             kmrnaDeg*x(1);         % mRNA degradation
             kProteinProd*x(1);     % protein production
             kProteinDeg*x(2)]);    % protein degradation

% Number of simulations;
NSims = 500;                        
% Result array
X_Array = zeros(2,probeNTimes,NSims);

% Plotting 
reverseStr='';
figure(1)
col=copper(NSims);

% Run multiple simulations
for isim = 1:NSims
    % Run simulation
    X_Array(:,:,isim) = SSA_direct_function(stoichiometry,prop,speciesCounts0,probeTimeArray);
    
    %%  <-- this gives a cell structure.  Hit %enter to run cell.
    
    % Plotting
    subplot(3,1,1)
    
    plot(probeTimeArray,X_Array(1,:,isim),'-*','color',col(isim,:));  % Plot mRNA vs time
    hold on
    
    subplot(3,1,2)
    plot(probeTimeArray,X_Array(2,:,isim),'-*','color',col(isim,:));  % Plot Protein vs time
    hold on
    
    subplot(3,1,3)
    plot(X_Array(1,:,isim),X_Array(2,:,isim),'-*','color',col(isim,:));  % Plot mRNA
    hold on
    
    % drawnow (print progress message to user)
    msg = sprintf('Processed %d/%d', isim, NSims);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end



%%
f2 = figure(2);
BINS = [linspace(0,30,40);linspace(0,200,40)];

for iT = 1:length(probeTimeArray)
f2 = figure(2); clf;
    for iSpec = 1:2
        subplot(2,1,iSpec)
        counts = hist(squeeze(X_Array(iSpec,iT,:)),BINS(iSpec,:));
        plot(BINS(iSpec,:),counts);
    end
    drawnow
    Movie(iT) = getframe(f2);
end

movie2avi(Movie,'tmp.avi')
    
    
