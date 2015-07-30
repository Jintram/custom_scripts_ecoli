%% stochastic version of nature model
%
% But I think it's easier to make a (very little) model myself; 
%
% Note furthermore the below code is horribly wrong, as it doesn't properly
% include small time steps. (I started writing without knowing where to go
% :)).



















%% stochastic version of nature model
REPEATS = 10;
TIMEPOINTS = 100;
rng(42); % seed

%% Parameters

% listed coefficients
intensityP      = 1;        % Does he mean E or p?
intensityMu     = .15;
intensityG      = .22;

transmissionEG      = 1.3;
transmissionMuE     = .7;
transmissionEE      = 1;

% Arbitrarily set (because only affects normalization)
% This is surely valid when the cross-correlation function is fitted
% directly. (See eq 13 and next paragraph suppl. Kiviet2014.)
% - For transmissionMuG, intensityMu sets the noise in traceMu, so also here
% valid.
% - For tranmissionEMu, ...?
transmissionMuG = 1;
tranmissionEMu = -1;


%% 

traceNoiseP  = rand(REPEATS,TIMEPOINTS) * intensityP;
traceNoiseMu = rand(REPEATS,TIMEPOINTS) * intensityMu;
traceNoiseG  = rand(REPEATS,TIMEPOINTS) * intensityG;

% note that everything is normalized in the paper.. 
% check how to appropriately deal with this

%% Now solve step by step
startE = 1;
startMu = .1;
startP = .1;
traceE=ones(REPEATS,TIMEPOINTS)*startE;
traceMu=ones(REPEATS,TIMEPOINTS)*startMu;
traceP=ones(REPEATS,TIMEPOINTS)*startP;
for R = 1
for t = 2:TIMEPOINTS
    
    % Instead of the model formula
    %traceE(R,t)  = traceP(R,t-1) - traceMu(R,t-1)*traceE(R,t-1);

    % A development is used, which looks like:
    traceE(R,t)  =                     traceP(R,t-1) + ... 
                     tranmissionEMu *  traceMu(R,t-1) -  ...
                                       traceE(R,t-1);
   % But I have no idea how this is derived, and how the normalization is 
   % done.

    traceMu(R,t) =   transmissionMuE * traceE(R,t-1) + ...
                     transmissionMuG * traceNoiseG(R,t-1) + ...
                                       traceNoiseMu(R,t-1); 

    traceP(R,t)  =    transmissionEE * traceE(R,t-1) + ... 
                      transmissionEG * traceNoiseG(R,t-1) + ...
                                       traceNoiseP(R,t-1);
       
    [traceE(R,t) , traceMu(R,t) ,  traceP(R,t)]
                                         
end
end
%% plot

figure, hold on;
%plot(traceE(1,1:25),'k')
plot(traceMu(1,1:25),'r')
plot(traceP(1,1:25),'b')

%% Calculate and plot some illustrative plots

% ===
figure, hold on;
plot(traceNoiseP(1,:), 'r-', 'Linewidth', 3);
xlabel('Time (a.u.)')
ylabel('Production (a.u.)')
FONTSIZE=20; set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal'); set(gca,'FontSize',FONTSIZE);

integratedTraceNoiseP1 = cumsum(traceNoiseP(1,:));

% ===
figure, hold on;
plot(integratedTraceNoiseP1, 'r-', 'Linewidth', 3);
xlabel('Time (a.u.)')
ylabel('Total # enzyme (a.u.)')
FONTSIZE=20; set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal'); set(gca,'FontSize',FONTSIZE);

smoothConcentration = smooth(integratedTraceNoiseP1,10);
plot(smoothConcentration, 'k--', 'Linewidth', 3);

% ===
figure, hold on;
cellularSize = smoothConcentration';
simpleApproximationTrace = integratedTraceNoiseP1./cellularSize;
plot(simpleApproximationTrace, 'b-', 'Linewidth', 3);
xlabel('Time (a.u.)')
ylabel('Enzyme / proxy cell size')
FONTSIZE=20; set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal'); set(gca,'FontSize',FONTSIZE);










