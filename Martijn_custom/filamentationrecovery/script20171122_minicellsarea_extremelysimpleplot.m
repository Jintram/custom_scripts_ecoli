

figure();

ylim([0,1]);
xlim([0,40]);

L=2:40;
miniRatio1=.5./L;
miniRatio2=(L-.5)./L;

hold on;
plot(L,miniRatio1)
plot(L,miniRatio2)