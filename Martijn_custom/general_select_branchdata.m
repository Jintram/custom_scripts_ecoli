function branchData = general_select_branchdata(branchData,starttrim,endtrim)
% function branchData = general_trim_branchdata(branchData,ntrim)
% 
% Cuts of first starttrim and last endtrim values from arrays for all 
% fields in branchData, throws away rest.

for idx = 1:length(branchData)
       
    % Trim all fields in current_branchData
    myfieldnames = fieldnames(branchData);
    for f_idx = 1:length(myfieldnames)
        current_field = myfieldnames(f_idx);
        branchData(idx).(current_field{1}) = branchData(idx).(current_field{1})(1+starttrim:end-endtrim);
    end
    
    % NaN detection (b/c trimming is usually done to remove NaN data).
    if any(isnan(branchData(idx).muP11_all)) || any(isnan(branchData(idx).time)) 
        disp('Still nan values detected!')
    end
    
end

end