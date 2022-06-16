%% Insert data

data=struct

% Intermediate [cAMP] (top)
data(1).name = 'top';
data(1).concentration = .314; %mM
data(1).s70 = [0.57, 0.53] % s70
data(1).CRP = [0.61, 0.7, 0.5, 0.68, 0.62] % CRP

% Higher [cAMP]
data(2).name = 'high';
data(2).concentration = 1.3; % mM
data(2).s70 = [0.82, 0.8, 0.73] %s70
data(2).CRP = [0.76, 0.82, 0.81, 0.77] %CRP

% Lower [cAMP]
data(3).name = 'low';
data(3).concentration = .05; %mM
data(3).s70 = [0.13, 0.4, 0.16] %s70
data(3).CRP = [0.18, 0.2] %CRP

%% And plot
figure(1); clf; hold on;

for i=1:numel(data)
        plot(...
            ones(1,numel(data(i).CRP))*data(i).concentration,...
            data(i).CRP,'ob','LineWidth',3,'MarkerSize',15)
        plot(...
            ones(1,numel(data(i).s70))*data(i).concentration,...
            data(i).s70,'s','LineWidth',3,'MarkerSize',15,'Color',[.4 .4 .4])
end

set(gca,'xscale','log');

title('Colony growth rates');
xlabel('cAMP [mM]');
ylabel('growth rate [dbl/hr]');
MW_makeplotlookbetter(20);


