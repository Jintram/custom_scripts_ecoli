
% Instructions
% ===
% 1. Run phospho_base_noise.m
% 2. Run 1st section of "phospho_analyze_branches.m" first to obtain
% branch structure.
%
% Structure of this file
% - First section analyzes one branch to illustrate effects of different
% weighing and normalizing of autocorr. (Using both matlab cross corr
% function, custom (MW) cross corr with weights, and DJK cross corr.)
% - Second section performs actual calcalation of autocorr, using
% DJK_getCrossVoc.m. 

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

%% Start of section II

%--------------------------------------------------------------------------
%                                                                         %
%        Section II - obtain autocorrelation for multiple strains         %
%                                                                         %
%-------------------------------------------------------------------------%

%% Now, get actual DJK composite autocorrs 

strain_names = fieldnames(myPhosphoData);

acs_DJK_Ryy={}; acs_DJK_t={};
for s_idx = 1:length(strain_names)

    current_strain = strain_names(s_idx);
    current_strain = current_strain{1};
    
    replicate_names = fieldnames(myPhosphoData.(current_strain))    
    for r_idx = 1:length(replicate_names)

        % Retrieve data from my data structure and prepare for DJK
        current_rep = replicate_names(r_idx);
        current_rep = current_rep{1};
        current_branchData = myPhosphoData.(current_strain).(current_rep).branchData;    
        current_p = myPhosphoData.(current_strain).(current_rep).p;
        current_p.timeField = 'time';

        % Trim data (remove last 1 values), b/c NaN values are present at end positions (?! TODO)
        current_branchData = general_trim_branchdata(current_branchData,1);    

        % Execute DJK cross corr
        [branches, crossCov_composite] = DJK_getCrossCov(current_p, current_branchData, 'muP11_all', 'muP11_all','extraNorm',1,'weighing',2)

        % Save output
        acs_DJK_t.(current_strain).(current_rep) = crossCov_composite.X;
        acs_DJK_Ryy.(current_strain).(current_rep) = crossCov_composite.Y./crossCov_composite.Y(1);

        % Tell user we're making progress
        disp(['Replicate ',current_rep,' done..']); 
    end
    disp(['Strain ',current_strain,' done..']); 
end
disp('All done');

currenttime = [date, '_', datestr(now, 'HH-MM-SS')];
save(['C:\Users\wehrens\Desktop\autocorr_',currenttime,'.mat']);

% TODO
% Additional trim for s734, r2 is actually required.
%{
current_strain='s734'; current_rep = 'r2';
current_branchData = myPhosphoData.(current_strain).(current_rep).branchData;
current_branchData = general_select_branchdata(current_branchData,1,1);
current_p = myPhosphoData.(current_strain).(current_rep).p;
current_p.timeField = 'time';
[branches, crossCov_composite] = DJK_getCrossCov(current_p, current_branchData, 'muP11_all', 'muP11_all','extraNorm',1,'weighing',2)
% Save output
acs_DJK_t.(current_strain).(current_rep) = crossCov_composite.X;
acs_DJK_Ryy.(current_strain).(current_rep) = crossCov_composite.Y./crossCov_composite.Y(1);
% Tell user we're making progress
disp(['Replicate ',current_rep,' done..']);
%}

%% And plot those
% WHAT DO YOU WANT TO PLOT??
% -----
current_strain='s735';
linecolorcounter=3;
% -----

replicate_names = fieldnames(myPhosphoData.(current_strain))

%% Preparing average plots

% Now to get error bars we have to perform some creative averaging
% - Bin points by timewindows
% - Base timewindows on highest dt

% Get some parameters necessary for the averaging
% - binwidth based on min dt (points should be independent!)
% - max tau on min time
dts=[]; end_times=[];
for r_idx = 1:length(replicate_names)
    % So get all timewindows    
    dts = [dts general_check_delta_t(acs_DJK_t.(current_strain).(current_rep))];
    current_times = acs_DJK_t.(current_strain).(current_rep);
    end_times = [end_times current_times(end)];
end
t_binwidth = min(dts)
max_tau = min(end_times)

Ryy_mean = []; Ryy_SEM = [];
bins_t_left_boundary = [t_binwidth/2:t_binwidth:max_tau-t_binwidth/2];
used_bin_centers = [];
datapoints_per_bin = [];
for t = bins_t_left_boundary

    data_in_bin = [];
    % Collect data 
    for r_idx = 1:length(replicate_names)
        
        current_rep = replicate_names(r_idx);
        current_rep = current_rep{1};
        
        % Find indices of data within bin
        indexes = find((acs_DJK_t.(current_strain).(current_rep)>t) .* (acs_DJK_t.(current_strain).(current_rep)<t+t_binwidth));
        
        current_Ryy = acs_DJK_Ryy.(current_strain).(current_rep);
        % Collect corresponding Ryy values to parameter data_in_bin
        data_in_bin = [data_in_bin, current_Ryy(indexes)];
    end
    
    % calculate means for this bin (only record when data available)
    if (length(data_in_bin)>0)
        % filter out nans
        data_in_bin = data_in_bin(find(~isnan(data_in_bin)))
        used_bin_centers(end+1) = (t+.5*t_binwidth);
        %data_in_bin        
        Ryy_mean(end+1) = mean(data_in_bin);
        Ryy_SEM(end+1) = std(data_in_bin)/sqrt(length(data_in_bin)); 
        datapoints_per_bin(end+1) = length(data_in_bin);
    end
end

% figure(1); plot(Ryy_mean,'o');
% figure(3); hold on;
% plot([Ryy_mean+Ryy_SEM],'.r'); 
% plot([Ryy_mean-Ryy_SEM],'.b'); 

% Get cell cycle times
mean_interDivTimes = [];
for r_idx = 1:length(replicate_names)
    
    current_rep = replicate_names(r_idx);
    current_rep = current_rep{1};
    
    interDivTimes = [myPhosphoData.(current_strain).(current_rep).s_all.interDivTime]; % {1} is to convert
    
    non_nans_idxs = find(~isnan(interDivTimes)); % get rid of NaNs
    mean_interDivTimes = [mean_interDivTimes, mean(interDivTimes(non_nans_idxs))];
end
mean_interDivTimes


%% Plot individual lines
% --------------------
mymaxtime = 300; % limit of x-axis
% --------------------

figure(1);
clf; hold on; legendtext={};

% also plot
t_maxes=[];
for r_idx = 1:length(replicate_names)
    current_rep = replicate_names(r_idx); % for title
    current_rep = current_rep{1};
    
    % plot
    plot(acs_DJK_t.(current_strain).(current_rep),acs_DJK_Ryy.(current_strain).(current_rep),['-'],'LineWidth',3,'Color',preferredcolors(r_idx,:));  %,somemarkers(r_idx)
    legendtext{end+1}=['ac rep ', current_rep]; % {1} to convert > str
        
    % Just for plot scaling
    t_maxes(end+1)=max(acs_DJK_t.(current_strain).(current_rep));
end

% zero line
plot([0, max(t_maxes)],[0,0],'k'); 

for r_idx = 1:length(replicate_names)
    % cell cycle time for colony
    plot(mean_interDivTimes(r_idx),0,'^','MarkerFaceColor',preferredcolors(r_idx,:),'MarkerEdgeColor','k','MarkerSize',12);
    
end

% cell cycle time averaged over all colonies
plot(mean(mean_interDivTimes),0,'^','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',12);

% cosmetics
set(gca,'FontSize',20);
%legend(legendtext,'Location','Best');
xlabel('Time (minutes)');
ylabel('Correlation (normalized)');

% Set x-limit to 3x cellcycle
%xlim([min(used_bin_centers), 4*mean(mean_interDivTimes)]);
% Set x-limit to half of total time
xlim([0, mymaxtime]);

%% Plot
figure(1); clf; hold on;
hE=errorbar(used_bin_centers,Ryy_mean,Ryy_SEM,'Color',preferredcolors(linecolorcounter,:));
set(hE,'Marker','none');

% cell cycle time for colony
plot(mean_interDivTimes,zeros(1,length(mean_interDivTimes)),'^','MarkerFaceColor',[.5,.5,.5],'MarkerEdgeColor',[.5,.5,.5],'MarkerSize',12);
% cell cycle time averaged over all colonies
plot(mean(mean_interDivTimes),0,'^','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',12);

% Set x-limit to 3x cellcycle
xlim([min(used_bin_centers), 4*mean(mean_interDivTimes)]);

% Cosmetics
set(gca,'FontSize',20);
title('Average autocorrelation');
xlabel('\tau (min)');
ylabel('Autocorrelation (normalized)');
plot([min(used_bin_centers), max(used_bin_centers)],[0,0],'k'); % zero line






