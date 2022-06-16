

% Some additional settings
FIGNUMBER=1;

figure(FIGNUMBER);
set(gca, 'xscale', 'log');
xlim([100,1000]);
ylim([0,1.4]);

MW_makeplotlookbetter(20);
