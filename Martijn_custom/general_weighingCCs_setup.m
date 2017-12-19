
%%

EXPORTFOLDER='\\storage01\data\AMOLF\users\wehrens\Latex3\Thesis\Chapter2_Methods\Figures\MatlabExport\';

associatedFieldNames={'t'    'x'    'y'};
FIELDPREFIX = '';


%% This sets up the correct fields 
% if you're working with real data to directly feed into applicable
% .._getCrossCor's subsection

%branches = branch_groups
%p.selectionName=name_rm_branch

fieldX = [FIELDPREFIX associatedFieldNames{1,2}]
fieldY = [FIELDPREFIX associatedFieldNames{1,3}]

p.timeField =associatedFieldNames{1}
p.movieName = 'test';

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

hInput=figure(2); clf; hold on;
plot(testbranches(IDX).(p.timeField), testbranches(IDX).(fieldX),'LineWidth',2)
plot(testbranches(IDX).(p.timeField), testbranches(IDX).(fieldY),'--','LineWidth',2)
legend({'signal 1','signal 2'});
ylim([-1,1])

xlabel('Time (a.u.)');
ylabel('Signal (a.u.)');

MW_makeplotlookbetter(14,[],[6,4]);

%% 
IDXs=[1,2];
hInput2=figure(101); clf; hold on;

for i=1:2
    subplot(3,1,(i)); hold on;
    plot(testbranches(IDXs(i)).(p.timeField), testbranches(IDXs(i)).(fieldX),'LineWidth',2)
    plot(testbranches(IDXs(i)).(p.timeField), testbranches(IDXs(i)).(fieldY),'--','LineWidth',2)
    
    xlim([0,100]);
    ylim([-1,1])
    %MW_makeplotlookbetter(14,[],[6,4]);
    MW_makeplotlookbetter(14,[],[6,6]);
end

xlabel('Time (a.u.)');
ylabel('Signal (a.u.)');

legend({'signal 1','signal 2'});

%% Creating CC with self

%{
testbranches(1).(fieldX)=testbranches(1).(fieldY)
testbranches(2).(fieldX)=testbranches(2).(fieldY)
branches=testbranches;
%testbranches(3).(fieldX)=testbranches(4).(fieldY)
%testbranches(3).(fieldX)=testbranches(4).(fieldY)
%}

%, , ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',ONSCREEN); 


%% Calculate the cross-correlation
p.weighing=4;
p.sameLength=1;
[branchesCorr, crossCor_composite] = MW_getCrossCor(p, branches, fieldX, fieldY)
%[branchesCorr, crossCor_composite] = DJK_getCrossCor(p, branches, fieldX, fieldY)

%%
hCC=figure(3); clf; hold on;

theColors=linspecer(4);

plot(branchesCorr(1).X,branchesCorr(1).Y,'LineWidth',13,'Color',theColors(1,:))
plot(branchesCorr(2).X,branchesCorr(2).Y,'LineWidth',11,'Color',theColors(2,:))
plot(branchesCorr(3).X,branchesCorr(3).Y,'LineWidth',9,'Color',theColors(3,:))
plot(branchesCorr(4).X,branchesCorr(4).Y,'LineWidth',7,'Color',theColors(4,:))
plot(crossCor_composite.X,crossCor_composite.Y,'ko-','LineWidth',3);

xlabel('? (a.u.)');
ylabel('Correlation');
MW_makeplotlookbetter(14,[],[6,4]);
%xlim([-30,30])

%%

hCC2=figure(102); clf; hold on;

theColors=linspecer(4);

for i=1:2
    subplot(3,1,(i)); hold on;
    plot(branchesCorr(i).X,branchesCorr(i).Y,'LineWidth',2,'Color',theColors(1,:))
    
    %xlim([0,100]);
    ylim([-1.5,1.5])
    MW_makeplotlookbetter(14,[],[6,6]);
end

%plot(branchesCorr(2).X,branchesCorr(2).Y,'LineWidth',11,'Color',theColors(2,:))
%plot(branchesCorr(3).X,branchesCorr(3).Y,'LineWidth',9,'Color',theColors(3,:))
%plot(branchesCorr(4).X,branchesCorr(4).Y,'LineWidth',7,'Color',theColors(4,:))
subplot(3,1,3); hold on;
plot(crossCor_composite.X,crossCor_composite.Y,'k-','LineWidth',2);
ylim([-1.5,1.5])

xlabel('Delay (a.u.)');
ylabel('Correlation');
MW_makeplotlookbetter(14,[],[6,6]);

%%
if exist('SAVEPLEASE','var')
    saveas(hInput2,[EXPORTFOLDER 'SVG_technical_CCs_exampleInput.svg']);
    saveas(hInput2,[EXPORTFOLDER 'TIF_technical_CCs_exampleInput.tif']);
    saveas(hInput2,[EXPORTFOLDER 'FIG_technical_CCs_exampleInput.fig']);

    saveas(hCC2,[EXPORTFOLDER 'SVG_technical_CCs_exampleCC.svg']);
    saveas(hCC2,[EXPORTFOLDER 'TIF_technical_CCs_exampleCC.tif']);
    saveas(hCC2,[EXPORTFOLDER 'FIG_technical_CCs_exampleCC.fig']);
end





