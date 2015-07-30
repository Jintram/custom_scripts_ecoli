function [the_bin_centers, the_means, the_stds, the_noises] = xy_mean_std_per_window(xdata, ydata, nr_windows,equal_datapoints_per_bin)
% function [the_bin_centers, the_means, the_stds, the_noises] = xy_mean_std_per_window(xdata, ydata, nr_windows,equal_datapoints_per_bin)


% There are two ways of dividing the data; an equal amount of datapoints
% per bin, or equal sizes of the bins.

if equal_datapoints_per_bin
    % equal nr datapoints: generate evenly spaced indexes for N bins
    window_indexes = round(linspace(0, length(xdata), nr_windows));
else
    % equal sizes of bins: divide x-data in evenly spaced windows
    bin_boundaries = linspace(0, max(xdata), nr_windows);
end

% loop over windows
the_bin_centers = [];
the_means = []; the_stds = []; the_noises = [];
for i = 2:nr_windows  
   
    % x values boundaries of window
    if equal_datapoints_per_bin
        left_boundary = xdata(window_indexes(i-1)+1);
        right_boundary = xdata(window_indexes(i));        
    else
        left_boundary = bin_boundaries(i-1);
        right_boundary = bin_boundaries(i);        
    end
    
    % obtain y values within that window
    y_indexes = find(xdata > left_boundary & xdata < right_boundary);
    window_y_data = ydata(y_indexes);
    
    % for convenience: calculate bin center
    the_bin_centers = [the_bin_centers (left_boundary+right_boundary)/2];
    
    % calculate statistics per window
    current_mean = mean(window_y_data);
    current_std = std(window_y_data);
    
    % and add these to vectors
    the_means = [the_means current_mean];
    the_stds = [the_stds current_std];
    the_noises = [the_noises current_std/current_mean];
       
end




