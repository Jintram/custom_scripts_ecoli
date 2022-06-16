
% 
% This script simply illustrates how the adder leads to shorter division
% times; and how this is reflected on normal and log scale.
% 
% Note that indeed as pointed out by the reviewer, the area under the L
% curve should be constant for constant mu and dL.


%% Understanding simple things..

cell1L0=2;
cell2L0=10;
deltaL = 2;

t=[0:.1:1.1];

cell1L = cell1L0 .* 2 .^ (.9*t);
cell2L = cell2L0 .* 2 .^ (.9*t);

treshold1 = cell1L0 + deltaL;
treshold2 = cell2L0 + deltaL;

figure; clf; hold on;
plot(t, cell1L,'LineWidth',2)
plot(t, cell2L,'LineWidth',2)

% tresholds
plot([min(t),max(t)],[cell1L0 cell1L0],':k','LineWidth',2)
plot([min(t),max(t)],[treshold1 treshold1],':k','LineWidth',2)
plot([min(t),max(t)],[cell2L0 cell2L0],':k','LineWidth',2)
plot([min(t),max(t)],[treshold2 treshold2],':k','LineWidth',2)

xlim([0,max(t)]);
ylim([0,3*cell2L0]);

xlabel('Time (hrs)');
ylabel('Length (um)');
MW_makeplotlookbetter(20);

%%
set(gca,'yscale','log');
ylim([0.01,3*cell2L0]);










