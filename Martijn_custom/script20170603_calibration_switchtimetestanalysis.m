function script20170603_calibration_switchtimetestanalysis()

%%
% the tube length here was approx. 107cm, and flow rate 8ul/min.

% Dataset 1
IMGPATHSTART='H:\EXPERIMENTAL_DATA_2017\2017_04_18_switchTimeTestLong_only\switch_time_test\long_tube_test_one\pos1\pos1-g-';
FRAMERANGE = [398:462,464:503];
SWITCHTIME = datenum(2017,04,20,15,17,20);

[t1,y1,normy1]=script20160603_z_signalanalysis(IMGPATHSTART,FRAMERANGE,SWITCHTIME);

% Dataset2
IMGPATHSTART='H:\EXPERIMENTAL_DATA_2017\2017_04_18_switchTimeTestLong_only\switch_time_test\long_tube_test_two\pos1\pos1-g-';
FRAMERANGE = [1:102];
SWITCHTIME = datenum(2017,04,21,09,30,26);

[t2,y2,normy2]=script20160603_z_signalanalysis(IMGPATHSTART,FRAMERANGE,SWITCHTIME);
normy2prime = 1-normy2;

%% Find the halfpoints! Bring them to me unspoiled.

halfpoint1Idx=find([0 normy1<0.5] & [normy1>0.5 0]);
halfpoint2Idx=find([0 normy2prime<0.5] & [normy2prime>0.5 0]);


%%
figure(1); clf; hold on;
plot(t1,normy1,'o-','LineWidth',2);
plot(t2,normy2prime,'o-','LineWidth',2);

ylabel('Normalized fluorescence');
xlabel('Time (mins)');
MW_makeplotlookbetter(20);

xlim([min([t1,t2]),max([t1,t2])]);

plot(t1(halfpoint1Idx),normy1(halfpoint1Idx),'ko','MarkerSize',20,'LineWidth',2)
plot(t2(halfpoint2Idx),normy2prime(halfpoint2Idx),'ko','MarkerSize',20,'LineWidth',2)

switchDelay=mean([t1(halfpoint1Idx) t2(halfpoint2Idx)])

%% PART II -- AFTER CHANGING THE TUBE LENGTH
% the tube length here was 15cm, and also flow rate 8ul/min.

IMGPATHSTART='H:\EXPERIMENTAL_DATA_2017\2017_04_18_switchTimeTestLong_only\switch_time_test\short_tube_two_switches\pos1\pos1-g-';
FRAMERANGE = [103:135];
SWITCHTIME = datenum(2017,04,21,12,25,25);

[t3,y3,normy3]=script20160603_z_signalanalysis(IMGPATHSTART,FRAMERANGE,SWITCHTIME);

IMGPATHSTART='H:\EXPERIMENTAL_DATA_2017\2017_04_18_switchTimeTestLong_only\switch_time_test\short_tube_two_switches\pos1\pos1-g-';
FRAMERANGE = [116:162]; % note there's some overlap here
SWITCHTIME = datenum(2017,04,21,12,49,05);

[t4,y4,normy4]=script20160603_z_signalanalysis(IMGPATHSTART,FRAMERANGE,SWITCHTIME);
normy4prime=1-normy4;

%% Find halfpoints
halfpoint3Idx=find([0 normy3<0.5] & [normy3>0.5 0]);
halfpoint4Idx=find([0 normy4prime<0.5] & [normy4prime>0.5 0]);



%%
figure(2); clf; hold on;
plot(t3,normy3,'o-','LineWidth',2);
plot(t4,normy4prime,'o-','LineWidth',2);

ylabel('Normalized fluorescence');
xlabel('Time (mins)');
MW_makeplotlookbetter(20);

xlim([min([t3,t4]),max([t3,t4])]);

plot(t3(halfpoint3Idx),normy3(halfpoint3Idx),'ko','MarkerSize',20,'LineWidth',2)
plot(t4(halfpoint4Idx),normy4prime(halfpoint4Idx),'ko','MarkerSize',20,'LineWidth',2)

switchDelayShortTube=mean([t3(halfpoint3Idx) t4(halfpoint4Idx)])





end





