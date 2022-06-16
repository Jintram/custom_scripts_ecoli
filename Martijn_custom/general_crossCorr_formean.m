
% Simple script exploring how to plot cross corr for mean trace
% I.e. make a simple cross corr function R(tau) for two parameters, as is
% done in MW_delayedScatter.
% Note that for the cross corrs we normally calculate, the colony-mean for
% that timepoint is already subtracted, therefor, "global" fluctuations are
% already partially filtered out.
%
% Note that 

% Plot a trace
figure(1), plot(output.branchavg.G6_mean_cycCor,'LineWidth',3), ylim([0,5000])

% Create mean-substracted parameters
fieldname1 = 'G6_mean_cycCor';
fieldname2 = 'muP5_fitNew_cycCor';

mean1 = mean(output.branchavg.(fieldname1))
noise1 = (output.branchavg.(fieldname1) - mean1)/mean1;

mean2 = mean(output.branchavg.(fieldname2));
noise2 = (output.branchavg.(fieldname2) - mean2)/mean2;

% Plot those
figure(2); clf; hold on;
plot(noise1,'b')
plot(noise2,'r')
        
% Calculate and plot cross correlation
figure(3); clf; hold on;
   
MAXLAGS = numel(noise1);

[Rtau, tau] = xcorr(...
        noise1, ...
        noise2, ...
        MAXLAGS,'coeff');

plot(tau,Rtau, '-', 'LineWidth', 3)
ylim([-1 1]);





