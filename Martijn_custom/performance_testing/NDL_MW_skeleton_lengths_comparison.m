
load([p.tracksDir p.movieName '-skeletonData.mat']);

figure(1); clf; hold on; title('MW vs NDL');
plot([allLengthsOfBacteriaInMicrons{:}],[allLengthsOfBacteriaInMicronsMW{:}],'.');
plot([0,max([allLengthsOfBacteriaInMicrons{:},allLengthsOfBacteriaInMicronsMW{:}])],[0,max([allLengthsOfBacteriaInMicrons{:},allLengthsOfBacteriaInMicronsMW{:}])]);
xlabel('NDL'); ylabel('MW');
MW_makeplotlookbetter(15);

figure(2); clf; hold on; title('MW vs. area');
plot([allPixelAreaOfBacterium{:}],[allLengthsOfBacteriaInMicronsMW{:}],'.');
xlabel('area'); ylabel('MW');
MW_makeplotlookbetter(15);

figure(3); clf; hold on; title('NDL vs. area');
plot([allPixelAreaOfBacterium{:}],[allLengthsOfBacteriaInMicrons{:}],'.');
xlabel('area'); ylabel('NDL');
MW_makeplotlookbetter(15);