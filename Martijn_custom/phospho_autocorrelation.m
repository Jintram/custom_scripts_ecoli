
% Instructions
% ===
% 1. Run phospho_base_noise.m
% 2. Run 1st section of "phospho_analyze_branches.m" first to obtain
% branch structure.

some_colors;

current_branchData = myPhosphoData.('s732').('r1').branchData;

%%% TODO, this is very weird, last muP11 is NaN
current_branchData = general_trim_branchdata(current_branchData,1);

current_p = myPhosphoData.('s732').('r1').p;
branch_nr = 200;

% So we want an autocorrelation, e.g. of the following data:
myfr = current_branchData(branch_nr).frame_nrs;
myti = current_branchData(branch_nr).time;
myY = current_branchData(branch_nr).muP11_all;
myLambda = current_branchData(branch_nr).count; % number of time datapoint is used for in branch structure

% Check whether delta t values are equal
dts = general_check_delta_t(myti);
dt=mean(dts); % note that time steps SHOULD be equal!!

% get manual ac without correction for branches
[ac,ac_raw, w]=general_get_autocorrelation(myY,0);
ac_w = general_get_autocorrelation(myY,myLambda);

% get matlab ac
[ac_ml,lags]=xcov(myY,'coeff');
ml_times = lags.*dt;

%% Plot my autocors
figure(1);
clf; hold on; legendtext={};
plot([0 myti],ac/ac(1),'-b','LineWidth',2); legendtext{end+1}='ac';
plot([0 myti],ac_raw/ac_raw(1),'--b','LineWidth',2); legendtext{end+1}='ac raw';
plot([0 myti],ac_w/ac_w(1),'--g','LineWidth',2); legendtext{end+1}='ac w';

plot(ml_times,ac_ml,'-r','LineWidth',1); legendtext{end+1}='ac ml';

legend(legendtext,'Location','Best');
xlabel('Time (minutes)');
ylabel('Correlation (normalized)');


% Get DJK autocorrelation for single branch, to compare
current_p.timeField = 'time';
singleBranch = current_branchData(branch_nr);
[branches, crossCov_composite] = DJK_getCrossCov(current_p, singleBranch, 'muP11_all', 'muP11_all','extraNorm',1,'weighing',2)
ac_DJK = crossCov_composite.Y./crossCov_composite.Y(1);

% Plot it too
plot(crossCov_composite.X,ac_DJK,'-k','LineWidth',2); legendtext{end+1}='ac DJK';
legend(legendtext,'Location','Best');
xlabel('Time (minutes)');
ylabel('Correlation (normalized)');

%% Now, get actual DJK composite autocorrs 

replicate_names = fieldnames(myPhosphoData.('s732'))
acs_DJK_Ryy={}; acs_DJK_t={};
for r_idx = 1:length(replicate_names)
        
    % Retrieve data from my data structure and prepare for DJK
    current_rep = replicate_names(r_idx);
    current_branchData = myPhosphoData.('s732').(current_rep{1}).branchData;    
    current_p = myPhosphoData.('s732').('r1').p;
    current_p.timeField = 'time';
    
    % Trim data (remove last 1 values), b/c NaN values are present at end positions (?! TODO)
    current_branchData = general_trim_branchdata(current_branchData,1);    
    
    % Execute DJK cross corr
    [branches, crossCov_composite] = DJK_getCrossCov(current_p, current_branchData, 'muP11_all', 'muP11_all','extraNorm',1,'weighing',2)
    
    % Save output
    acs_DJK_t{end+1} = crossCov_composite.X;
    acs_DJK_Ryy{end+1} = crossCov_composite.Y./crossCov_composite.Y(1);
    
    % Tell user we're making progress
    disp(['Replicate ',current_rep{1},' done..']); % {1} to convert > str
end

disp('All done');

%% And plot those

figure(1);
clf; hold on; legendtext={};

% also plot
for r_idx = 1:length(replicate_names)
    current_rep = replicate_names(r_idx); % for title
    
    % plot
    plot(acs_DJK_t{r_idx},acs_DJK_Ryy{r_idx},'-xk','LineWidth',1); 
    legendtext{end+1}=['ac rep ', current_rep{1}]; % {1} to convert > str
        
end

% cosmetics
legend(legendtext,'Location','Best');
xlabel('Time (minutes)');
ylabel('Correlation (normalized)');


%% 

% Now to get error bars we have to perform some creative averaging
% - Bin points by timewindows
% - Base timewindows on highest dt

% Get some parameters necessary for the averaging
% - binwidth based on min dt (points should be independent!)
% - max tau on min time
dts=[]; end_times=[];
for r_idx = 1:length(replicate_names)
    % So get all timewindows    
    dts = [dts general_check_delta_t(acs_DJK_t{r_idx})];
    current_times = acs_DJK_t{r_idx};
    end_times = [end_times current_times(end)];
end
t_binwidth = min(dts)
max_tau = min(end_times)

Ryy_mean = []; Ryy_SEM = [];
bins_t_left_boundary = [t_binwidth/2:t_binwidth:max_tau-t_binwidth/2];
bins_t_center = bins_t_left_boundary+t_binwidth/2
for t = bins_t_left_boundary

    data_in_bin = [];
    % Collect data 
    for r_idx = 1:length(replicate_names)
        % Find indices of data within bin
        indexes = find((acs_DJK_t{r_idx}>t) .* (acs_DJK_t{r_idx}<t+t_binwidth));
        
        current_Ryy = acs_DJK_Ryy{r_idx};
        % Collect corresponding Ryy values to parameter data_in_bin
        data_in_bin = [data_in_bin, current_Ryy(indexes)];
    end
    
    % calculate means for this bin
    
    Ryy_mean(end+1) = mean(data_in_bin);
    Ryy_SEM(end+1) = std(data_in_bin)/sqrt(length(data_in_bin)); 
end

% figure(1); plot(Ryy_mean,'o');
% figure(3); hold on;
% plot([Ryy_mean+Ryy_SEM],'.r'); 
% plot([Ryy_mean-Ryy_SEM],'.b'); 

% Get cell cycle times
mean_interDivTimes = [];
for r_idx = 1:length(replicate_names)
    
    current_rep = replicate_names(r_idx);
    
    interDivTimes = [myPhosphoData.('s732').(current_rep{1}).s_all.interDivTime]; % {1} is to convert
    
    non_nans_idxs = find(~isnan(interDivTimes)); % get rid of NaNs
    mean_interDivTimes = [mean_interDivTimes, mean(interDivTimes(non_nans_idxs))];
end
mean_interDivTimes

%% Plot
figure(1); clf; hold on;
hE=errorbar(bins_t_center,Ryy_mean,Ryy_SEM,'Color',preferredcolors(1,:));
set(hE,'Marker','none');

% cell cycle time for colony
plot(mean_interDivTimes,zeros(1,length(mean_interDivTimes)),'^','MarkerFaceColor',[.5,.5,.5],'MarkerEdgeColor',[.5,.5,.5],'MarkerSize',12);
% cell cycle time averaged over all colonies
plot(mean(mean_interDivTimes),0,'^','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',12);

% Set x-limit to 3x cellcycle
xlim([min(bins_t_center), 4*mean(mean_interDivTimes)]);

% Cosmetics
set(gca,'FontSize',20);
title('Average autocorrelation');
xlabel('\tau (min)');
ylabel('Autocorrelation (normalized)');
plot([min(bins_t_center), max(bins_t_center)],[0,0],'k'); % zero line






