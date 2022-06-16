% plots several qualitatively cell cycle behaviours (especially for
% concentration)
% CASE (5) IS YFP of MG22

t=0:0.001:1; % phase
V=2.^t;     % Volume
meanVol=mean(V);   % average Volume in case you want to normalize
figure(1)   % Volume & total protein
figure(2)   % Protein concentration
figure(3)   % Production rate

% (1) same #gene copies through whole cell cycle (hypothetic)
r1=ones(1,1001); % rate
P1=1+t;   % total protein
% (2) gene doubling at phase=0.2
r2=zeros(1,1001);
r2(1:201)=0.5/0.9;
r2(202:1001)=1/0.9;
P2=ones(1,1001);
for i=1:1001;
P2(i)=sum(r2(1:i-1))*0.001+1;
end
% (3) gene doubling at phase=0.5
r3=zeros(1,1001);
r3(1:501)=0.5/0.75;
r3(502:1001)=1/0.75;
P3=ones(1,1001);
for i=1:1001;
P3(i)=sum(r3(1:i-1))*0.001+1;
end
% (4) gene doubling at phase=0.8
r4=zeros(1,1001);
r4(1:801)=0.5/0.6;
r4(802:1001)=1/0.6;
P4=ones(1,1001);
for i=1:1001;
P4(i)=sum(r4(1:i-1))*0.001+1;
end
% (5) cst rate till phase=0.5 then linear increase until doubling at
% phase=1 (CFP/YFP) (version 2013-01-07)
r5=zeros(1,1001);
r5(1:501)=0.8;
r5(501:1001)=0.8+0.8*((t(501:1001)-t(501))/0.5);
P5=ones(1,1001);
for i=1:1001;
P5(i)=sum(r5(1:i-1))*0.001+1;
end
% % (5) cst rate till phase=***0.6*** then linear increase until doubling at
% % phase=1 (CFP/YFP)
%r5=zeros(1,1001);
%r5(1:601)=1/1.2;
%r5(602:1001)=1/1.2+1/1.2*((t(602:1001)-t(602))/0.4);
%P5=ones(1,1001);
%for i=1:1001;
%P5(i)=sum(r5(1:i-1))*0.001+1;
%end
% (6) linear increasing rate till phase=0.6 (->rate doubled), then cst (GFP
% of ASC631 2012-06-17)
r6=zeros(1,1001);
r6(1:601)=1/1.7+1/1.7*(t(1:601)/0.6);
r6(602:1001)=2/1.7;
P6=ones(1,1001);
for i=1:1001;
P6(i)=sum(r6(1:i-1))*0.001+1;
end
% (7) linear increase with phase until doubling. typical phase correction
% (Rosenfeld2005)
r7=zeros(1,1001);
r7=2/3+2/3*t;
P7=ones(1,1001);
for i=1:1001;
P7(i)=sum(r7(1:i-1))*0.001+1;
end

figure(1) %***** TOTAL LENGTH AND FLUO *****
% if divided by 'meanVol' -> average volume =1
clf
%plot(t,V/meanVol,'k','LineWidth',2);
plot(t,V/meanVol,'k','LineWidth',2);
hold on
title('total')
%plot(t,P1,'Color', [1 0.8 0.2],'LineWidth',2)  % cst copy #
%plot(t,P2,'Color',[1 0 0],'LineWidth',2)       % doubling 0.2
%plot(t,P3,'Color',[0.8 0 0.8],'LineWidth',2)    % doubling 0.5
%plot(t,P4,'Color',[0.4 0 0.8],'LineWidth',2)  % doubling 0.8
plot(t,P5/meanVol,'Color',[0 0.7 1],'LineWidth',2)     % YFP/CFP: cst, then lin increase
%plot(t,P6,'Color',[0 0.7 0])     % GFP: lin increase, then cst
%plot(t,P7,'Color',[0.7 0 0.7])   % lin increase (std cecy corr)

%ylim([0.95 2.05])
%legend('Vol','cst copy','0.2','0.5','0.8','Location','NW')
legend('Vol','Location','NW')
%legend('Vol','YFP/CFP','GFP','lin increase','Location','NW')
%legend('Vol','YFP/CFP','Location','NW')

figure(2)  %****CONCENTRATION****
clf
hold on
title('conc')
%plot(t,P1./V,'r')
%plot(t,P2./V,'Color',[1 0 0],'LineWidth',2)       % doubling 0.2
%plot(t,P3./V,'Color',[0.8 0 0.8],'LineWidth',2)    % doubling 0.5
%plot(t,P4./V,'Color',[0.4 0 0.8],'LineWidth',2)  % doubling 0.8
%plot(t,P1./V,'Color', [1 0.8 0.2],'LineWidth',2)
%plot(t,P2./V,'Color',[1 0 0],'LineWidth',2)
%plot(t,P3./V,'Color',[0.8 0 0.8],'LineWidth',2)
%plot(t,P4./V,'Color',[0.4 0 0.8],'LineWidth',2)
plot(t,P5./V,'Color',[0 0.7 1],'LineWidth',2)  
%plot(t,P6./V,'Color',[0 0.7 0])
%plot(t,P7./V,'Color',[0.7 0 0.7])

figure(3)  %****RATE*****
%if divided by initial value, all init prod rates are set equal to =1.
clf
hold on
title('rate')
%plot(t,r1,'Color', [1 0.8 0.2],'LineWidth',2)
%plot(t,r2/r2(1),'Color',[1 0 0],'LineWidth',2)
%plot(t,r3/r3(1),'Color',[0.8 0 0.8],'LineWidth',2)
%plot(t,r4/r4(1),'Color',[0.4 0 0.8],'LineWidth',2)
plot(t,r5,'Color',[0 0.7 1],'LineWidth',2)  
%plot(t,r6,'Color',[0 0.7 0])
%plot(t,r7,'Color',[0.7 0 0.7])


