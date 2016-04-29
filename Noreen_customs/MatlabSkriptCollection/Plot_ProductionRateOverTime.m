% plots average production rate over time and standard deviation

% *** adjust ******
pstruct=p;  %required for name of figure!
myschnitzcells=s_rm_fitTime;
timefield='time_atdR';
ratefield='dR5_cycCor';
% *****************

ratetotal=[];
timetotal=[];

for i=1:length(myschnitzcells)
    s=myschnitzcells(i);
    if s.useForPlot==1
        if length(s.(timefield))>0 & length(s.(timefield))==length(s.(ratefield))
            ratetotal=[ratetotal; s.(ratefield)'];
            timetotal=[timetotal; s.(timefield)'];
        end
    end
end

% get average per time point. cells weighed equally.
timeunique=unique(timetotal);
meanrate=[];
stdevrate=[];
timeTolerance=0.5; % sometimes images taken a little later/earlier
for i=1:length(timeunique)
    time=timeunique(i);
    idxt=find(timetotal>time-timeTolerance & timetotal<time+timeTolerance);
    meanrate_at_t=mean(ratetotal(idxt));
    stdrate_at_t=std(ratetotal(idxt));
    meanrate=[meanrate, meanrate_at_t];
    stdevrate=[stdevrate, stdrate_at_t];
end

figure
clf
title([pstruct.movieDate ' ' pstruct.movieName ' . ' ratefield], 'Interpreter','None');
hold on
errorbar(timeunique,meanrate,stdevrate)
xlabel('time (min)')
ylabel([ratefield], 'Interpreter','None')