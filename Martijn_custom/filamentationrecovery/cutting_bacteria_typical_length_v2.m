% MW 01/07/2014




%% Some bacterial length for which to generate data
lengths=[3:0.1:50];
% Typical length at which division takes place (due to waves of min)
typical = 3.77/2;

% admin
l_out = []; cuts_at = []; cuts_at_fraction = [];
% loop over lengths
for l=lengths
       
    % calculate number of cuts possible with this typical length
    nrcuts = round(l/typical);
    lresult = l/nrcuts;
    
    % represent this data as Rutger did
    for i=1:nrcuts-1 
        % just x-axis for plot
        l_out = [l_out, l];
        % calculate where cut i takes place for this length l
        current_cut_at = lresult*i;
        % add to vector (not used/plotted)
        cuts_at   = [cuts_at, current_cut_at];
        % calculate fraction at which cut is, to plot on y-axis
        cuts_at_fraction = [cuts_at_fraction, current_cut_at/l];
    end
    
end

%% plot
figure(3); clf;

%%  load data Rutger
load('D:\Local_Software\Martijn_extensions\Martijn_custom\filamentationrecovery\RR_asymmetry_pattern_all.mat')

hold on;
plot(A(:,1)',A(:,2)','.','Color',[.4 .4 .4])

%% plot
plot(l_out, cuts_at_fraction,'o');
title('Dividing bacteria');
xlabel('Length');
%ylabel('Fraction at which division takes place');
ylabel('L_d/L_m');
MW_makeplotlookbetter(20);





