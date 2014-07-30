
% Instructions
% ===
% 1. Run phospho_base_noise.m
% 2. Run 1st section of "phospho_analyze_branches.m" first to obtain
% branch structure.


current_branchData = myPhosphoData.('s732').('r1').branchData;

%%% TODO, this is very weird, last muP11 is NaN
ntrim = 1; disp('Trimming data.')
for idx = 1:length(current_branchData)
       
    % Trim all fields in current_branchData
    myfieldnames = fieldnames(current_branchData);
    for f_idx = 1:length(myfieldnames)
        current_field = myfieldnames(f_idx);
        current_branchData(idx).(current_field{1}) = current_branchData(idx).(current_field{1})(1:end-ntrim);
    end
    
    if any(isnan(current_branchData(idx).muP11_all)) || any(isnan(current_branchData(idx).time)) 
        disp('Still nan values detected!')
    end
    
end

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


%% Get DJK autocorrelation
current_p.timeField = 'time';
current_branchData = current_branchData(branch_nr); % TODO REMOVE - select only 1
[branches, crossCov_composite] = DJK_getCrossCov(current_p, current_branchData, 'muP11_all', 'muP11_all','extraNorm',1,'weighing',2)
ac_DJK = crossCov_composite.Y./crossCov_composite.Y(1);

%% Plot it.
%figure(2); clf; hold on;
plot(crossCov_composite.X,ac_DJK,'-k','LineWidth',2); legendtext{end+1}='ac DJK';
legend(legendtext,'Location','Best');
xlabel('Time (minutes)');
ylabel('Correlation (normalized)');







