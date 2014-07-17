% phoshpo_analyze_branches

branchData = DJK_getBranches(p,s_all,'dataFields',{'time','muP11_all'});

figure; clf; hold on;
plot(branchData(1).time, branchData(1).muP11_all);

sliding_length_vs_noise