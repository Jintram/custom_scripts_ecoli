


% divide everything by biomass
% (aka everything is expressed as fraction biomass)
%
%        Ec-gene |<------
%         |             |
%        \ /            |
%        Ec             |      Ribo = 1-Ec
% sugar --------> metabolite ------> biomass
%           |           |
%           -------------
%           fast inhibition (post transcript.)

INITIALPROTIENVALUE     = 0;
AIMEDPROTEINFRACTION    = .1;   % = Keq
LAMBDA                  = 2;    % /hr; ln(2)*mu=.7*mu=lambda
NORMALIZEDCELLMASS0     = 1;

prodrate = LAMBDA;

dt=.01;


proteinTrace = [INITIALPROTIENVALUE];
for t = [1:dt:10]
    
    deltaProtein        = dt * (prodrate - LAMBDA)
    0                   = dt * (prodrate - AIMEDPROTEINFRACTION LAMBDA)
    proteinTrace(end+1) = proteinTrace(end) + deltaProtein;
    
    
end



figure(1); clf; hold on;
plot(proteinTrace,'.');

ylim([0,1]);

xlabel('Time');
ylabel('C-sector as fraction of total protein mass');






