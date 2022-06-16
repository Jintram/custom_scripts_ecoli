
%function [p,schnitzcells] = MW_plotting_mu_t_raw(multiple_ps) 
% Function that returns a scatter plot of all mu(t) values, with mu on
% y-axis and t on x-axis. 

% Gather information
alltimes = ...
    [schnitzcells.time];

allmus = ...
    [schnitzcells.muP11_fitNew_all];

% Scatterplot
figure(1);
clf();
plot(alltimes,allmus,'-x')
xlabel('t')
ylabel('\mu')