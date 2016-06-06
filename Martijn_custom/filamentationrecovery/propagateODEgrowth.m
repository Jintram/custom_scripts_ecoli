
%%
parameters.mu = 1;
parameters.divisionTime = 1; % for timer model to work, this needs to be equal to mu
parameters.divisionSize = 6;
parameters.addedSize = 6;
parameters.divisionType = 'timer';

% noise parameters
parameters.fluctuationIntensity = 0.04; % 0.04 seems good value
parameters.relaxationTimeFluctuationsMinutes = 50; % 50 seems good value

INITIALLENGTH=25;

times = 1:400; % simulation time
divisionEvents = [times(1)]; divisionSizes=[INITIALLENGTH];
growthEfficiencies = [1];
%%
cellLengths=[INITIALLENGTH];
for time = times(2:end)
    [cellLengths,divisionEvents,divisionSizes,growthEfficiencies] = ODEgrowth(cellLengths,time,divisionEvents,divisionSizes,growthEfficiencies,parameters);
end

%%

figure(1); clf; hold on;
plot(times,cellLengths,'LineWidth',3)

plot(divisionEvents,zeros(1,numel(divisionEvents)),'^r','MarkerFaceColor','r','MarkerSize',15)

xlabel('time (min)');
ylabel('cellength [\mum]');
MW_makeplotlookbetter(20);

figure(2); clf; hold on;
plot(times,growthEfficiencies,'LineWidth',3)

plot(divisionEvents,zeros(1,numel(divisionEvents)),'^r','MarkerFaceColor','r','MarkerSize',15)

xlabel('time (min)');
ylabel('growth efficiency');
MW_makeplotlookbetter(20);
