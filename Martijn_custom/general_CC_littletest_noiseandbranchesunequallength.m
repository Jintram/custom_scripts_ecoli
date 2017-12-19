%%
%
% Branches should be normalized before cross-correlations are determined.
% This is done subtracting the colony mean at that point in time. This is
% illustrated by the below script.
%
% Note that the cross-correlation script had to be modified to be able to
% handle branches of unequal length, but this was not the case for the
% normalization script.
%
% - Martijn Wehrens, 2017-11
%

clear testbranches2

lineColors=linspecer(3);

N=30;
Nprime=15;

t=[1:N];
idxsPrime=[1:Nprime]+7;

trend= 2*t*1/50;

testbranches2(1).count=ones(1,N)*2;
testbranches2(1).schnitzNrs=ones(1,N);
testbranches2(1).t=t;
testbranches2(1).X=rand(1,N)+trend;
testbranches2(1).Y=rand(1,N);

testbranches2(2).count=ones(1,N)*2;
testbranches2(2).schnitzNrs=ones(1,N);
testbranches2(2).t=t;
testbranches2(2).X=rand(1,N)+trend;
testbranches2(2).Y=rand(1,N);

p.dataFields={'t','X','Y'}

%%

testbranches2 = DJK_addToBranches_noise(p, testbranches2)


%%
h1=figure(1); clf; 
subplot(2,2,1); hold on;
plot(testbranches2(1).t,testbranches2(1).X,'Color',lineColors(1,:));
plot(testbranches2(2).t,testbranches2(2).X,'Color',lineColors(2,:));
xlabel('Time (a.u.)');
ylabel('Signal (a.u.)');

subplot(2,2,2); hold on;
plot(testbranches2(1).t,testbranches2(1).noise_X,'Color',lineColors(1,:));
plot(testbranches2(2).t,testbranches2(2).noise_X,'Color',lineColors(2,:));
xlabel('Time (a.u.)');
ylabel('Signal (a.u.)');


%% Now add a branch with unequal length and repeat

testbranches2(3).count=ones(1,Nprime)*2;
testbranches2(3).schnitzNrs=ones(1,Nprime);
testbranches2(3).t=t(idxsPrime);
testbranches2(3).X=rand(1,Nprime)+trend(idxsPrime);
testbranches2(3).Y=rand(1,Nprime);

testbranches2 = DJK_addToBranches_noise(p, testbranches2)

%%

subplot(2,2,3); hold on;
plot(testbranches2(1).t,testbranches2(1).X,'-','Color',lineColors(1,:));
plot(testbranches2(2).t,testbranches2(2).X,'-','Color',lineColors(2,:));
plot(testbranches2(3).t,testbranches2(3).X,'Color',lineColors(3,:),'LineWidth',2);
xlabel('Time (a.u.)');
ylabel('Signal (a.u.)');

subplot(2,2,4); hold on;
plot(testbranches2(1).t,testbranches2(1).noise_X,'-','Color',lineColors(1,:));
plot(testbranches2(2).t,testbranches2(2).noise_X,'-','Color',lineColors(2,:));
plot(testbranches2(3).t,testbranches2(3).noise_X,'Color',lineColors(3,:),'LineWidth',2);
xlabel('Time (a.u.)');
ylabel('Signal (a.u.)');

MW_makeplotlookbetter(14,[],[12.8,8]);




