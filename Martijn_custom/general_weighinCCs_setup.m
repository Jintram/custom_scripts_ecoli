
%% This sets up the correct fields 
% if you're working with real data to directly feed into applicable
% .._getCrossCor's subsection

%branches = branch_groups
%p.selectionName=name_rm_branch

fieldX = [FIELDPREFIX associatedFieldNames{1,2}]
fieldY = [FIELDPREFIX associatedFieldNames{1,3}]

p.timeField =associatedFieldNames{1}

%or...
%fieldX='noise_C6_mean_cycCor'; fieldY='muP9_fitNew_cycCor'; p.timeField = 'C_time'

%% This will provide you with some branches to test the cross-correlation function

N=100;

testbranches=struct;

testbranches.schnitzNrs = ones(1,N);
testbranches.(p.timeField)   = 1:N;
testbranches.(fieldX) = ...
    repmat([ones(1,10)*-.5 ones(1,10)*.5],[1,N/(2*10)]);
testbranches.(fieldY) = ...
    repmat([ones(1,5)*-.5 ones(1,10)*.5 ones(1,5)*-.5],[1,N/(2*10)]);
testbranches.count = ones(1,N);
testbranches.branchpoints= zeros(1,N);

branches=testbranches;

%% And add a second branch only 3/5 as long

N=40;
testbranches(2).schnitzNrs = ones(1,N);
testbranches(2).(p.timeField)   = 1:N;
testbranches(2).(fieldX) = ...
    repmat([ones(1,10)*-.5 ones(1,10)*.5],[1,N/(2*10)]);
testbranches(2).(fieldY) = ...
    repmat([ones(1,5)*-.5 ones(1,10)*.5 ones(1,5)*-.5],[1,N/(2*10)]);
testbranches(2).count = ones(1,N);
testbranches(2).branchpoints= zeros(1,N);

branches=testbranches;

%% And add a third and fourth branch that show redundancy in the counts

for idx=3:4
    N=100; % needs to be dividable by 20
    testbranches(idx).schnitzNrs = ones(1,N);
    testbranches(idx).(p.timeField)   = 1:N;
    testbranches(idx).(fieldX) = ...
        repmat([ones(1,10)*-.5 ones(1,10)*.5],[1,N/(2*10)]);
    testbranches(idx).(fieldY) = ...
        repmat([ones(1,5)*-.5 ones(1,10)*.5 ones(1,5)*-.5],[1,N/(2*10)]);
    testbranches(idx).count = [2.*ones(1,N/2) ones(1,N/2)];
    
    testbranches(idx).branchpoints = [0.*ones(1,N/2) ones(1,N/2)];
end

branches=testbranches;

%%
IDX=1;

figure(2); clf; hold on;
plot(testbranches(IDX).(p.timeField), testbranches(IDX).(fieldX),'LineWidth',2)
plot(testbranches(IDX).(p.timeField), testbranches(IDX).(fieldY),'--','LineWidth',2)
legend({'signal 1','signal 2'});
ylim([-1,1])

xlabel('Time');
ylabel('Signal');

MW_makeplotlookbetter(20);

%% Creating CC with self

testbranches(1).(fieldX)=testbranches(1).(fieldY)
testbranches(2).(fieldX)=testbranches(2).(fieldY)
branches=testbranches;
%testbranches(3).(fieldX)=testbranches(4).(fieldY)
%testbranches(3).(fieldX)=testbranches(4).(fieldY)








%, , ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',ONSCREEN); 