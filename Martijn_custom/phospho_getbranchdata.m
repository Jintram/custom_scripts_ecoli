function myPhosphoData = phospho_getbranchdata(myPhosphoData,strain,rep)

    disp(['Calculating branch data for ',strain,', rep ',rep,'..']);

    p = myPhosphoData.(strain).(rep).p;
    s_all = myPhosphoData.(strain).(rep).s_all;

    if isfield(s_all,'mu11_fitNew_all')
        myPhosphoData.(strain).(rep).theMuField = 'mu11_fitNew_all';
    elseif isfield(s_all,'muP11_all')
        myPhosphoData.(strain).(rep).theMuField = 'muP11_all'
    else
        disp('This is a problem, myPhosphoData, error #1');
    end
    
    % note that old data has field 'frames' instead of 'frame_nrs'
    if isfield(s_all,'frames')
        myPhosphoData.(strain).(rep).branchData = DJK_getBranches(p,s_all,'dataFields',{'time',myPhosphoData.(strain).(rep).theMuField,'length_fitNew','frames'});
        for idx = 1:length(myPhosphoData.(strain).(rep).branchData)
            myPhosphoData.(strain).(rep).branchData(idx).frame_nrs = myPhosphoData.(strain).(rep).branchData.frames;
        end
    else
        myPhosphoData.(strain).(rep).branchData = DJK_getBranches(p,s_all,'dataFields',{'time',myPhosphoData.(strain).(rep).theMuField,'length_fitNew','frame_nrs'});
    end
    
    %{
    % After we have data available in separate branches, we want to
    % calculate auto correlation.
    
    branches = DJK_addToBranches_noise(p, branchData,'dataFields',{'dR5_time'  'R_time'  'muP11_fitNew_atdR5' 'muP11_fitNew_atdR5_cycCor' 'dR5_cycCor'  'dG5_cycCor' 'dG5' 'dR5' });
    trimmed_branches = DJK_trim_branch_data(branches,4);
    branch_groups = DJK_divide_branch_data(trimmed_branches);

    DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dG5_cycCor', 'noise_muP11_fitNew_atdR5_cycCor','selectionName',name_rm_branch,'timeField','R_time');
    %}

end