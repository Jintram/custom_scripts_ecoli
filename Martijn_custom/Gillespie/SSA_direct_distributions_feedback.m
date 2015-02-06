
% Written by Johannes Keegstra, adapted by Martijn Wehrens
% 2015/02

some_colors;
Nspecies = 5;

% Specify initial and final times
time    = 0;
tstop   = 60; % in seconds

% Initial concentrations of species -- COULD be initialized with s.s.
% (TODO? cq to test, reaching s.s. is also good test?)
speciesCounts0      = zeros(Nspecies,1); 

% Stoichiometry matrix (e.g. 2A -> B)
stoichiometry = ...
... %   rxn1,   rxn2,  rxn3,  rxn4,  rxn5,  rxn6,  rxn7,  rxn8,  rxn9,  rxn10, rxn11    
    [  [ 1  ,  -1  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,  -1  ,   0  ,   0  ];...    % a   #1
       [ 0  ,   1  ,  -1  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,  -1  ,   0  ];...    % b   #2 
       [ 0  ,   0  ,   1  ,  -1  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,  -1  ];...    % c   #3
       [ 0  ,   0  ,   0  ,   0  ,   1  ,  -1  ,   0  ,   0  ,   0  ,   0  ,   0  ];...    % Eab #4
       [ 0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   1  ,  -1  ,   0  ,   0  ,   0  ] ...    % Ebc #5
    ];
    % Reactions: 0 -> a; a -> b; b -> c; c -> 0; 0 -> Exx; Exx -> 0.

% Define times at which simulation is halted to probe concentrations 
probeNTimes = 100;
probeTimeArray = linspace(0,tstop,probeNTimes);  % times at which to plot.

% General constants
muHr          = .93; % doubling rate / hr - Ref:Jain2009
lambdaHr      = log(2)*muHr; % exponential growth rate / hr
lambdaSec     = lambdaHr/3600; % lambda [/s]
    % Aimed steady state E values
    Ess     = 3000; % enzyme steady state level, based on 3mln * 1000 ppm.
    mss     = 100; % metabolite steady state level; this number is interesting as it will determine (kab and thus) sensitivity perturbations
% Reaction constants
kImport     = 562066/2 % [/s] importing nutrient - glucose (see ppt), and Neidhardt1996; assume 1/2 since 1:1 growth/ATP
kab         = kImport/(Ess*mss) % rxn a>b; assuming full flux /s.; ignore dilution OK?
kbc         = kImport/(Ess*mss) % rxn b>C
kDilution   = lambdaSec % dilution
kConsumeC   = kImport/mss-lambdaSec % consumption rxn C
kProdEab    = kDilution/Ess % production Eab
kProdEbc    = kDilution/Ess % production Ebc

% Calculate propensities
prop = @(x)([kImport;            % rxn1, import gluc
             kab * x(4) * x(1);  % rxn2, convert a->b, catal. Eab       
             kbc * x(5) * x(2);  % rxn3, convert b->c, catal. Ebc
             kConsumeC  * x(3);  % rxn4, consume c for growth
             
             kProdEab;           % rxn5, make enzyme
             kDilution * x(4);   % rxn6, dilute enzyme
             kProdEbc;           % rxn7, make enzyme
             kDilution * x(5);   % rxn8, dilute enzyme
             
             kDilution * x(1);   % rxn9, dilute a
             kDilution * x(2);   % rxn10, dilute b
             kDilution * x(3);   % rxn11, dilute c
             ]);    

% Number of simulations;
NSims = 3;                        
% Result array
X_Array = zeros(2,probeNTimes,NSims);

% For message
reverseStr='';
% Plotting 
figure(1)
col=copper(NSims);

% Run multiple simulations
for isim = 1:NSims
    disp(['Starting sim ' num2str(isim)]);
    
    % Run simulation
    X_Array(:,:,isim) = SSA_direct_function(stoichiometry,prop,speciesCounts0,probeTimeArray,Nspecies);
    
    %%  <-- this gives a cell structure.  Hit %enter to run cell.
    
    % Plotting
    subplot(2,1,1)
    
    plot(probeTimeArray,X_Array(1,:,isim),'-r','color',preferredcolors(2,:));  % Plot a
    plot(probeTimeArray,X_Array(2,:,isim),'-g','color',preferredcolors(3,:));  % Plot b
    plot(probeTimeArray,X_Array(3,:,isim),'-b','color',preferredcolors(4,:));  % Plot b
    hold on
    
    subplot(2,1,2)
    plot(probeTimeArray,X_Array(4,:,isim),'-r','color',preferredcolors(2,:));  % Plot a
    plot(probeTimeArray,X_Array(5,:,isim),'-g','color',preferredcolors(3,:));  % Plot b
    hold on       
    
    % plot(X_Array(1,:,isim),X_Array(2,:,isim),'-','color',col(isim,:));  % space of 1 vs. 2
    
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
    
    
