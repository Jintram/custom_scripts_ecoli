

function [means,stds] sliding_window_stds(tvalues,yvalues,nrbins,plotting);
% This function takes t values, y values as input (t because usually we're 
% looking at a time series), and determines for time windows of length
% max(tvalues)/nrbins the means and standard deviations (stds), and returns
% those.
%
% If plotting=1 then it will also plot the distributions it obtains.


nr_bins = 10;
mymap = colorGray(nr_bins);
yvalue_means = [];
yvalue_stds = [];

% figure
figure(1);
clf;
hold on;

% timewindows
maxt=max(plotx_732_pos1);
dt=maxt/nr_bins;
mytwindows=[0:dt:maxt];

for idx_window = 1:(length(mytwindows)-1)
  
    disp(num2str(idx_window));
    
    % get data for this window
    value_indices = find(plotx_732_pos1>mytwindows(idx_window) & plotx_732_pos1 < mytwindows(idx_window+1));
    yvalues_window = ploty_732_pos1(value_indices);
    
    size(yvalues_window)
    
    % perform binning again (make fn for that!)
    [yscores_window] = hist(yvalues_window,mybins);
    
    % normalizing - TODO needed/what about weighing ?!?!
    A732=sum(yscores_window)*dx;
    yscores_window=yscores_window./A732;
    
    % plotting
    plot(mybins, yscores_window,'-','LineWidth',2,'color',mymap(idx_window,:));%,'color',colorAmolfGreen); 
    
    % get statistics
    yvalue_means(end+1) = mean(yvalues_window);
    yvalue_stds(end+1) = std(yvalues_window);
    
end

yvalue_means
yvalue_stds

figure(2)
plot(yvalue_stds)
ylim([0,max(yvalue_stds)*1.1])


end