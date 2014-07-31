% NW_trim_branch_data_time trims data from beginning of branches untill the
% start time of the branches is larger a given initial time
%
function trimmed_branches = NW_trim_branch_data_time(branches,timefield,inittime);

% GET BRANCH INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% investiage first branch and find idx of the time points larger than
% inittime
timevec=branches(1).(timefield);
idx=find(timevec>inittime);

% TRIM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trimmed_branches = branches;
fields = fieldnames(branches);
for f = 1:length(fields)
  for i=1:length(trimmed_branches)
    trimmed_branches(i).(char(fields(f))) = branches(i).(char(fields(f)))(idx);
  end
end