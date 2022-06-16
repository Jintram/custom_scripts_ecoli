% plots sum of total fluor over time

timefield='G_time';
%timefield='G_frames';
fluofield='R5_sum';
fluoaveragefield='G5_mean';
lengthfield='length_fitNew';

alltime=[schnitzcells.(timefield)];
allfluo=[schnitzcells.(fluofield)];
allavfluo=[schnitzcells.(fluoaveragefield)];
alllength=[schnitzcells.(lengthfield)];

timeunique=unique(alltime);
sumfluo=zeros(size(timeunique));
errorfluo=zeros(size(timeunique));
meanfluo=zeros(size(timeunique));
sumlength=zeros(size(timeunique));

for i=1:length(timeunique)
    idx=find(alltime==timeunique(i));
    sumfluo(i)=sum(allfluo(idx));
    errorfluo(i)=std(allfluo(idx));
    meanfluo(i)=mean(allavfluo(idx));
    sumlength(i)=sum(alllength(idx));
end

figure(4)
set(gcf,'WindowStyle','docked')
errorbar(timeunique,sumfluo,errorfluo)
title ('all schnitz data, no weighing')

%figure
%plot(timeunique,meanfluo.*sumlength)
%title('conc * length')