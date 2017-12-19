


%% 
% Note that Wallden et al. assume that the duration of the C+D period
% depends on the growth rate.
%
% Importantly, Ho & Amir (2015) present a model for initiation of 
% replication which has the following components
% - the cell produces initiators + autorepressors
% - initiatiors are broken down after repl. initiation, which means that
%   the number of initiators is equal to the added volume to the cell
% - for next initiation, an additional # of initiators is needed that
%   saturates all oris', hence delta=ori*delta0

% I typically used perturbations of 2X up/down, or for larger ones 20X
% up/down. For ori's I used factor of 20X up, 8X down.

SAVEDIR = 'D:\Local_Data\Dropbox\Dropbox\Filamentation recovery\Submission_Current_Biology_round2\RebuttalFigs\AmirModel\';
if ~exist('scriptExecCounter')
    scriptExecCounter = 0;
end

%%
% simulation parameters
simpleVersion=1;
silence=1; % suppress messages to user
perturbationTime = 25; % note the code under that if-statement should be adjusted to your needs!

dt = 0.001;
tau=80/60; % I used 20 or 80 typically
mu = 1/tau; % 1/tau=mu; mu in dbl/hr
S0=2; % initiation volume
ENDTIME=250;
CPERIODTIME=40/60; % Cooper & Helmstetter
DPERIODTIME=20/60;
deltaS=2;
%deltaS = 1/(2^((CPERIODTIME+DPERIODTIME)/tau));

% initialization
CperiodCounters=[];
DperiodCounters=[];

% initial conditions
t=0;
chromos=1;
%S=20*2 
%oris=20*2^floor((CPERIODTIME+DPERIODTIME)/tau);
%{
oris=2*2^floor((CPERIODTIME+DPERIODTIME)/tau);
S=deltaS*oris;
CperiodCounters=[];
for idx=1:floor(((CPERIODTIME+DPERIODTIME)/tau))
    CperiodCounters = [idx*tau CperiodCounters];
end
Sreference=S;
%}

% Initial conditions from eqiulibration run
% For 3dbl/hr
%{
S = 21.7897;
oris =    16;
CperiodCounters = [0.8150    0.4810    0.1470];
SreferenceEnd = 16.0893;
%}
% For 1dbl/hr
S= 2.7345;
oris=2;
CperiodCounters =  0.2680;
Sreference = 2.3802;


% more initialization
counter=0;

%% Calculate some expected quantities
expectedSize = deltaS*2^((CPERIODTIME+DPERIODTIME)*mu);
disp(['Expected size=' num2str(expectedSize)]);
expectedOris = 2^((CPERIODTIME+DPERIODTIME)*mu); 
disp(['Expected ori''s=' num2str(expectedOris)]);
% chromo expectation = 1, not counting stuff that is not fully replicated


%%
disp('Simulation starting.');
disp('----------------');
%debugtracker=0;
perturbationPerformedAlready = 0;
while t<ENDTIME
   
    %%
    % counter
    counter=counter+1;
    
    % update time
    t(end+1)        = t(end)        + dt;
    % update cell size
    S(end+1)        = S(end)        + S(end)*log(2)*mu*dt;
    % keep chromo #, but might be updated later
    chromos(end+1)  = chromos(end);
    % keep ori #, but might be updated later
    oris(end+1)       = oris(end);
      
    % update all ongoing replication cycles
    for CDidx = 1:numel(CperiodCounters)
        CperiodCounters(CDidx) = CperiodCounters(CDidx) + dt;                
    end
    
    % initiate a replication process if required and double oris
    if (S(end)-Sreference) > (deltaS*oris(end))
        
        % initiate new replication counter
        CperiodCounters(end+1)=dt; 
        
        % double the number of ori's
        oris(end)=oris(end)*2;
        
        % re-start making 'initiator molecule', i.e. reference size
        Sreference=S(end);
        
        % user message
        if ~silence
            disp(['Replication initiated at t=' num2str(t(end)) ', current delta is ' num2str(deltaS*oris(end)) '']);
        end
        %disp(['>> '];
        
    end
    
    if simpleVersion
        
        % if replication/division cycle is done, double the amount of chromosomes,
        % remove timer, and start division process
        toRemove=[];
        for CDidx = 1:numel(CperiodCounters)
            if CperiodCounters(CDidx)>(CPERIODTIME+DPERIODTIME)
                
                if ~silence
                    disp(['Cell divided at t=' num2str(t(end)) ', delta is now = ' num2str(deltaS*oris(end))]);
                end
                
                % Actual division process
                S(end) = S(end)/2; % divide cell size
                oris(end) = oris(end)/2; % divide ori's                
                
                % Update reference size to memorize also the lost volume
                Sreference=Sreference-S(end);
                
                % double chromo's
                % chromos(end) = chromos(end)*2; 
                
                % this timer is now done, note to remove it later
                toRemove(end+1)=CDidx; 
            end
        end
        CperiodCounters(toRemove) = []; % remove timer
        
    else
        
        % if replication cycle is done, double the amount of chromosomes,
        % remove timer, and start division process
        toRemove=[];
        for CDidx = 1:numel(CperiodCounters)
            if CperiodCounters(CDidx)>CPERIODTIME
                chromos(end) = chromos(end)*2; % double chromo's            
                DperiodCounters(end+1)=dt; % initiate division delay
                toRemove(end+1)=CDidx; % note which timer to remove later            
                if ~silence
                    disp(['Replication ended, division started at t=' num2str(t(end))]);
                end
            end
        end
        CperiodCounters(toRemove) = []; % remove timer

        % update all ongoing divisions 
        for CDidx = 1:numel(DperiodCounters)
            DperiodCounters(CDidx) = DperiodCounters(CDidx) + dt;        
        end


        % if division cycle is done, divide S, chromos and oris by 2
        % remove timer, and start division process
        toRemove=[];
        for CDidx = 1:numel(DperiodCounters)
            if DperiodCounters(CDidx)>DPERIODTIME
                % Actual division process
                chromos(end) = chromos(end)/2; % divide chromo's
                S(end) = S(end)/2; % divide cell size
                oris(end) = oris(end)/2; % divide ori's
                toRemove(end+1) = CDidx; % note which timer to remove later
                if ~silence
                    disp(['Cell divided at t=' num2str(t(end)) ', delta is now = ' num2str(deltaS*oris(end))]);
                end

                % 
                Sreference=S(end);
            end
        end
        DperiodCounters(toRemove) = []; % remove timer
        
    end    
    
    if t(end)>perturbationTime & ~perturbationPerformedAlready
        S(end)=S(end)*8;
        oris(end)=oris(end)*8;
        %oris(end)=oris(end)/8;
        perturbationPerformedAlready=1;
        disp(['Perturbing the cell at t=' num2str(t(end))]);
    end
    
    %debugtracker(end+1) = S(end)-Sreference;
end
disp('----------------');
disp('Simulation done');

%% 
disp('Summary of end state');
disp('--------------');
SEnd=S(end)
tEnd=t(end)
orisEnd=oris(end)
CperiodCountersEnd=CperiodCounters
SreferenceEnd=Sreference
disp('--------------');


%%


h1=figure(1); clf; hold on;
% data
plot(t,S,'LineWidth',2)
plot(t,oris,'-','LineWidth',2)
% expectations
plot([min(t), max(t)],[expectedSize expectedSize],'--','LineWidth',2);
plot([min(t), max(t)],[expectedOris expectedOris],'--','LineWidth',2);
%plot(t,chromos,':','LineWidth',2)
xlabel('Time (hrs)');
ylabel('Size, S (um)');

ylim([0, max([S,oris]*1.1)]);
xlim([min(t), max(t)]);

%plot(t,debugtracker,'-');

%legend({'Expected size','Expected ori''s','Size','Chromosomes','Ori''s'});
%legend({'Expected size','Expected ori''s','Size','Ori''s'});
legend({'Size','Ori''s','Expected size','Expected ori''s'});

title(['Mu = ' num2str(mu) ' dbl/hr']);

MW_makeplotlookbetter(20);


%%
disp('Saving figures');
SIZE=[5.5,5.5];
OFFSET = [2,2];
set(h1,'Units','centimeters','Position',[OFFSET SIZE]*2)
MW_makeplotlookbetter(8*2);%,optionalParameters);

set(h1,'RendererMode','manual','Renderer','Painters');

scriptExecCounter = scriptExecCounter +1; % just for naming params.

saveas(h1, [SAVEDIR 'export_' num2str(scriptExecCounter) '.svg']);
saveas(h1, [SAVEDIR 'export_' num2str(scriptExecCounter) '.fig']);

%%
disp('All done');










