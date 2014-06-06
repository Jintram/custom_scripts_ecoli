
function [t_center,means,stds,N] = sliding_window_stds(nrybins,tvalues,yvalues,nrtbins,plotting,fignumstart,mylinecolor);
% This function takes t values, y values as input.
% It divides the data up in timewindows of width max(tvalues)/nrbins.
% For each of these windows, the y values are binned according to ybins
% (these are bins with width dy),
% and the means and standard deviations (stds) of these distributions are
% also determined.
%
% If plotting=1 then it will also plot the distributions it obtains.
% t_center returns the centers of the time windows used.

totaldy = max(yvalues)-min(yvalues);
dy=totaldy/nrybins;
ybins=[min(yvalues):dy:max(yvalues)];

%{
mymap = [[1, 0, 0];[0, 1, 0];[0, 0, 1]];
if nrtbins>3 % fn only seems to work for >2 bins
    mymap = colorGray(nrtbins);
end
%}
t_center = [];
means = [];
stds = [];
N = [];

% figure
if plotting
    figure(fignumstart);
    hold on;
    set(gca,'FontSize',20);
    title('Distribution of doubling time');
    xlabel('\mu (doublings/hr)');
    ylabel('Probability (normalized)');
end

% timewindows
maxt=max(tvalues);
dt=maxt/nrtbins;
mytwindows=[0:dt:maxt];

for idx_window = 1:(length(mytwindows)-1)
     
    % get data for this window
    value_indices = find(tvalues>mytwindows(idx_window) & tvalues < mytwindows(idx_window+1));
    yvalues_window = yvalues(value_indices);    
    
    % perform binning again (make fn for that!)
    [yscores_window] = hist(yvalues_window,ybins);
    
    % normalizing - TODO needed/what about weighing ?!?!
    area=sum(yscores_window)*dy;
    yscores_window=yscores_window./area;
    
    % plotting
    if plotting
        plot(ybins, yscores_window,'-','LineWidth',2,'color',mylinecolor);%,'color',colorAmolfGreen); 
    end
    
    % get statistics
    t_center(end+1) = (mytwindows(idx_window)+mytwindows(idx_window+1))/2;
    means(end+1) = mean(yvalues_window);
    stds(end+1) = std(yvalues_window);
    N(end+1) = length(yvalues_window);
    
end

end