% IMPORTANT:
% CHANGE THE #BINS, PERIODICITY, NOISE SUBTRACTION (CONC) BELOW WITHIN THE
% SCRIPT

%% RATES YFP CFP
ph=[];
yfp=[];
cfp=[];
time=[];
schn=[];
frames=[];
frames_unique=[];
highlightschnitzes=[]; %[100:40:400 210 182 423 ];%[363 446 398]; %blubb
% TAKE: [ 210 182 423]  (2012-05-08)
%modifications: "blubb" -> remove these lines for "standard plotting"


myschnitzcells=s_rm_fitTime;


for i=1:length(myschnitzcells)
    s=myschnitzcells(i);
    if s.useForPlot==1
        if length(s.time_atdY)>0 & length(s.phase2_atdY)==length(s.dY5_sum_dt) & length(s.phase2_atdY)==length(s.time_atdY) &  s.completeCycle==1
            ph=[ph; s.phase2_atdY']; %BLUBB!!! wrong phase/time maybe! !!! fits to dY5_sum_dt
            yfp=[yfp; s.dY5'];%[yfp; s.dY5_sum_dt'];%[yfp; s.dY5_sum_dt_s'];
            cfp=[cfp; s.dC5'];%[cfp; s.dC5_sum_dt'];%[yfp; s.dY5_sum_dt_s'];
            %yfp=[yfp; s.dY5'];
            %cfp=[cfp; s.dC5'];
            %yfp=[yfp; s.dY5_sum_dt'];
            %cfp=[cfp; s.dC5_sum_dt'];
            %yfp=[yfp; s.dY5_sum_dt'/mean(s.dY5_sum_dt)];
            %cfp=[cfp; s.dC5_sum_dt'/mean(s.dC5_sum_dt)];
            reltime=s.time_atdY-s.time(1);
            time=[time; reltime'];
            schnrep=zeros(1,length(s.phase2_atdY))+i;
            schn=[schn; schnrep'];
            framesrep=zeros(1,length(s.phase2_atdY))+length(s.frames);
            frames=[frames; framesrep'];
            frames_unique=[frames_unique; length(s.frames)];
        end
    end
end

% normalize and save unnormalized vector
yfpabs=yfp;
cfpabs=cfp;
yfp=yfp/mean(yfp);
cfp=cfp/mean(cfp);

% binning RATES YFP CFP
%idx=find(frames<61 & frames>54);
idx=find(frames<665 & frames>0); % loose conditions like frames >0 and frames<200 include all cells
%plot(ph(idx),yfp(idx),'.r')
%plot(ph(idx),cfp(idx),'.','Color',[0 0.7 0])
schnsub=schn(idx);
phsub=ph(idx);
yfpsub=yfp(idx);
cfpsub=cfp(idx);


% do binning
numbins=8;
binborders=[0:1/numbins:0.98];
binborders=[binborders,1];
binnedcfp=[];
binnedyfp=[];
meanbin=[0.5/numbins:1/numbins:0.98];
% stddev vector
binnedstdcfp=[];
binnedstdyfp=[];
% std error vector  (stddev/sqrt(#datapoints). don't use. see below
binnederrorcfp=[];
binnederroryfp=[];
% bootstrapping vector. use. see below.
binnedbootmeancfp=[];
binnedbootmeanyfp=[];

for i=1:numbins
    idx2=find(phsub>=binborders(i) & phsub<binborders(i+1));
    yfpinbin=yfpsub(idx2);
    cfpinbin=cfpsub(idx2);
    binnedyfp=[binnedyfp,mean(yfpinbin)];
    binnedcfp=[binnedcfp,mean(cfpinbin)];
    binnedstdyfp=[binnedstdyfp,std(yfpinbin)];
    binnedstdcfp=[binnedstdcfp,std(cfpinbin)];
    %NW 2013-12: I don't know if the following standard errors of the mean
    %may be used since we only talk about ONE experiment and not repeated
    %experiments:
    binnederroryfp=[binnederroryfp,std(yfpinbin)/sqrt(length(yfpinbin))];
    binnederrorcfp=[binnederrorcfp,std(cfpinbin)/sqrt(length(cfpinbin))];
    % NW2013-12: I believe that bootstrapping is applicable since it
    % resamples the distribution of the mean from one dataset. Condition is
    % that the individual data points are independent, which is likely
    % fullfilled because typically each cell has only one data point per
    % bin (subsequent datapoints might be dependent). 
    % Note: The resulting estimate of the error of mean is actually very similar
    % (quantitatively) than the SEM above (is this generally true??).
    currentbootmeancfp=bootstrp(10000,@mean,cfpinbin); %sample size 10,000. sample the 'mean'. use 'cfpinbin' as input distribution
    currentbootmeanyfp=bootstrp(10000,@mean,yfpinbin);
    binnedbootmeancfp=[binnedbootmeancfp, std(currentbootmeancfp)];  % standard deviation of the sampled mean
    binnedbootmeanyfp=[binnedbootmeanyfp, std(currentbootmeanyfp)];    
    clear cfpinbin yfpinbin currentbootmeancfp currentbootmeanyfp
end

figure(9)
clf
plot(ph(idx),yfp(idx),'.','Color',[1 0.7 0])%[0.9 0.8 0])
hold on
plot(meanbin,binnedyfp,'.-k','MarkerSize',15)
%plot(meanbin,4*binnedyfp,'.-','MarkerSize',10,'Color',[0.5 0.5 0.5])
%plot(meanbin,6*binnedyfp,'.-','MarkerSize',10,'Color',[0.5 0.5 0.5])
if ~isempty(highlightschnitzes)
    for j=1:length(highlightschnitzes)
       ss=highlightschnitzes(j);
       idxschn=find(schnsub==ss);
       phschn=phsub(idxschn);
       yfpschn=yfpsub(idxschn);
       plot(phschn,yfpschn,'.-','MarkerSize',15,'Color',[0.5 0.4 0])
    end
end
   
figure(10)
clf
plot(ph(idx),cfp(idx),'.','Color',[0 0.7 1])
hold on
plot(meanbin,binnedcfp,'.-k','MarkerSize',15)
%plot(meanbin,4*binnedcfp,'.-','MarkerSize',10,'Color',[0.5 0.5 0.5])
%plot(meanbin,6*binnedcfp,'.-','MarkerSize',10,'Color',[0.5 0.5 0.5])

figure(11)
clf
errorbar(meanbin,binnedyfp,binnedbootmeanyfp,'.-','MarkerSize',15,'Color',[1 0.7 0])%[0.9 0.8 0])
hold on
errorbar(meanbin,binnedcfp,binnedbootmeancfp,'.-','MarkerSize',15,'Color', [0 0.7 1])
grid on

% ********************************************************************************************
% ********************************************************************************************

%% CONC YFP CFP

phconc=[];
yfpconc=[];
cfpconc=[];
timeconc=[];
schnconc=[];
framesconc=[];
frames_uniqueconc=[];
abstimeconc=[];

myschnitzcells=s_rm_fitTime;
%highlight_schnitzes: hard-coded below

for i=1:length(myschnitzcells)
    s=myschnitzcells(i);
    if s.useForPlot==1
        if length(s.Y_time)>0 & length(s.phase2_atY)==length(s.Y5_mean) & length(s.phase2_atY)==length(s.Y_time) &  s.completeCycle==1
            phconc=[phconc; s.phase2_atY'];
            %yfpconc=[yfpconc; s.Y5_mean'/mean(s.Y5_mean)];
            %cfpconc=[cfpconc; s.C5_mean'/mean(s.C5_mean)];
            abstimeconc=[abstimeconc; s.Y_time'];
            yfpconc=[yfpconc; s.Y5_mean'];
            cfpconc=[cfpconc; s.C5_mean'];
            
            reltimeconc=s.Y_time-s.time(1);
            timeconc=[timeconc; reltime'];
            schnrepconc=zeros(1,length(s.phase2_atY))+i;
            schnconc=[schnconc; schnrepconc'];
            framesrepconc=zeros(1,length(s.phase2_atY))+length(s.frames);
            framesconc=[framesconc; framesrepconc'];
            frames_uniqueconc=[frames_uniqueconc; length(s.frames)];
        end
    end
end

% normalize and save unnormalized vector
yfpconcabs=yfpconc;
cfpconcabs=cfpconc;
yfpconc=yfpconc/mean(yfpconc);
cfpconc=cfpconc/mean(cfpconc);


% get noise values (population mean at each time point subtracted)
% all cells weighed equally
timeunique=unique(abstimeconc);
noiseyfpconc=zeros(size(yfpconc));
noisecfpconc=zeros(size(cfpconc));
meanyfp_vec=[];
meancfp_vec=[];
timeTolerance=0.5; % sometimes images taken a little later/earlier
for i=1:length(timeunique)
    time=timeunique(i);
    idxt=find(abstimeconc>time-timeTolerance & abstimeconc<time+timeTolerance);
    meanyfpconc_at_t=mean(yfpconc(idxt));
    meancfpconc_at_t=mean(cfpconc(idxt));
    noiseyfpconc(idxt)=yfpconc(idxt)-meanyfpconc_at_t;
    noisecfpconc(idxt)=cfpconc(idxt)-meancfpconc_at_t;
    meanyfp_vec=[meanyfp_vec, meanyfpconc_at_t];
    meancfp_vec=[meancfp_vec, meancfpconc_at_t];
end


figure(15)
clf
plot(timeunique,meanyfp_vec,'.-','MarkerSize',10,'Color',[1 0.7 0])%[0.9 0.8 0]) %don't divide by mean, at least not in 'noise'case
hold on
plot(timeunique,meancfp_vec,'.-','Color',[0 0.7 1],'MarkerSize',10)
grid on

% binning conc YFP CFP

% ***************  ADJUST ***************************
% choose if abs value (normalized) or noise (mean subtracted) is used
%yfpconcuse=yfpconc;  %BLUBB!!
%cfpconcuse=cfpconc;
yfpconcuse=noiseyfpconc;
cfpconcuse=noisecfpconc;
% ***************************************************

%idx=find(frames<61 & frames>54); % choose framerange
idx=find(frames<685 & frames>0); % loose conditions like frames >0 and frames<200 include all cells
schnsub=schn(idx);
phconcsub=phconc(idx);
cfpconcsub=cfpconcuse(idx);
yfpconcsub=yfpconcuse(idx);

% ********** ADJUST **************************
% choose if binning is supposed to be periodic (conc is continuous!)
periodic=0; % =1; =0; STANDARD=0
% *******************************************

if periodic
    
    % do binning periodic (pt at =0 and =1 are identical)
    numbins=6; %actually only (n-1) independent bins
    binborders=[0.5/(numbins-1):1/(numbins-1):0.98];
    meanbinconc=[0:1/(numbins-1):0.98];
    meanbinconc=[meanbinconc, 1];
    binnedyfpconc=[];
    binnedcfpconc=[];
    % stddev vector conc
    binnedstdcfpconc=[];
    binnedstdyfpconc=[];
    % std error vector conc (stddev/sqrt(#datapoints). don't use. see above.
    binnederrorcfpconc=[];
    binnederroryfpconc=[];
    % bootstrapping of mean. use. see above.
    binnedbootmeancfpconc=[];
    binnedbootmeanyfpconc=[];
    

    for i=1:(numbins-1)
        if i==1 % first=last idx, since periodic binning
            idx2=find(phconcsub<binborders(i) | phconcsub>=binborders(end));
            yfpconcinbin=yfpconcsub(idx2);
            cfpconcinbin=cfpconcsub(idx2);
            binnedyfpconc=[binnedyfpconc,mean(yfpconcinbin)];
            binnedcfpconc=[binnedcfpconc,mean(cfpconcinbin)];
            binnedstdyfpconc=[binnedstdyfpconc,std(yfpconcinbin)];
            binnedstdcfpconc=[binnedstdcfpconc,std(cfpconcinbin)];
            % standard error of mean. don't use. see above.
            binnederroryfpconc=[binnederroryfpconc,std(yfpconcinbin)/sqrt(length(yfpconcinbin))];
            binnederrorcfpconc=[binnederrorcfpconc,std(cfpconcinbin)/sqrt(length(cfpconcinbin))];
            % bootstrapping. use. see above. (description for prod. rates)
            currentbootmeancfpconc=bootstrp(10000,@mean,cfpconcinbin); %sample size 10,000. sample the 'mean'. use 'cfpconcinbin' as input distribution
            currentbootmeanyfpconc=bootstrp(10000,@mean,yfpconcinbin);
            binnedbootmeancfpconc=[binnedbootmeancfpconc, std(currentbootmeancfpconc)];  % standard deviation of the sampled mean
            binnedbootmeanyfpconc=[binnedbootmeanyfpconc, std(currentbootmeanyfpconc)];    
            clear cfpconcinbin yfpconcinbin currentbootmeancfpconc currentbootmeanyfpconc
        else
            idx2=find(phconcsub>=binborders(i-1) & phconcsub<binborders(i));
            yfpconcinbin=yfpconcsub(idx2);
            cfpconcinbin=cfpconcsub(idx2);
            binnedyfpconc=[binnedyfpconc,mean(yfpconcinbin)];
            binnedcfpconc=[binnedcfpconc,mean(cfpconcinbin)];
            binnedstdyfpconc=[binnedstdyfpconc,std(yfpconcinbin)];
            binnedstdcfpconc=[binnedstdcfpconc,std(cfpconcinbin)];
            % standard error of mean. don't use. see above
            binnederroryfpconc=[binnederroryfpconc,std(yfpconcinbin)/sqrt(length(yfpconcinbin))];
            binnederrorcfpconc=[binnederrorcfpconc,std(cfpconcinbin)/sqrt(length(cfpconcinbin))];
            % bootstrapping. use. see above. (description for prod. rates)
            currentbootmeancfpconc=bootstrp(10000,@mean,cfpconcinbin); %sample size 10,000. sample the 'mean'. use 'cfpconcinbin' as input distribution
            currentbootmeanyfpconc=bootstrp(10000,@mean,yfpconcinbin);
            binnedbootmeancfpconc=[binnedbootmeancfpconc, std(currentbootmeancfpconc)];  % standard deviation of the sampled mean
            binnedbootmeanyfpconc=[binnedbootmeanyfpconc, std(currentbootmeanyfpconc)];    
            clear cfpconcinbin yfpconcinbin currentbootmeancfpconc currentbootmeanyfpconc
        end
    end
    % at last binned data (periodic)
     binnedyfpconc=[binnedyfpconc,binnedyfpconc(1)];
     binnedcfpconc=[binnedcfpconc,binnedcfpconc(1)];
     binnedstdyfpconc=[binnedstdyfpconc,binnedstdyfpconc(1)];
     binnedstdcfpconc=[binnedstdcfpconc,binnedstdcfpconc(1)];
     binnederroryfpconc=[binnederroryfpconc,binnederroryfpconc(1)];
     binnederrorcfpconc=[binnederrorcfpconc,binnederrorcfpconc(1)];
     binnedbootmeancfpconc=[binnedbootmeancfpconc,binnedbootmeancfpconc(1)];
     binnedbootmeanyfpconc=[binnedbootmeanyfpconc,binnedbootmeanyfpconc(1)];

else   %STANDARD CASE
    % do binning not=periodic
    numbins=8;
    binborders=[0:1/numbins:0.98];
    binborders=[binborders,1];
    binnedyfpconc=[];
    binnedcfpconc=[];
    meanbinconc=[0.5/numbins:1/numbins:0.98];
    
    % stddev vector conc
    binnedstdcfpconc=[];
    binnedstdyfpconc=[];
    % std error vector conc (stddev/sqrt(#datapoints)  . don't use. see above.
    binnederrorcfpconc=[];
    binnederroryfpconc=[];
    % bootstrapping of mean. use. see above.
    binnedbootmeancfpconc=[];
    binnedbootmeanyfpconc=[];

    for i=1:numbins
        idx2=find(phconcsub>=binborders(i) & phconcsub<binborders(i+1));
        yfpconcinbin=yfpconcsub(idx2);
        cfpconcinbin=cfpconcsub(idx2);
        binnedyfpconc=[binnedyfpconc,mean(yfpconcinbin)];
        binnedcfpconc=[binnedcfpconc,mean(cfpconcinbin)];
        binnedstdyfpconc=[binnedstdyfpconc,std(yfpconcinbin)];
        binnedstdcfpconc=[binnedstdcfpconc,std(cfpconcinbin)];
        binnederroryfpconc=[binnederroryfpconc,std(yfpconcinbin)/sqrt(length(yfpconcinbin))];
        binnederrorcfpconc=[binnederrorcfpconc,std(cfpconcinbin)/sqrt(length(cfpconcinbin))];
        % bootstrapping. use. see above. (description for prod. rates)
        currentbootmeancfpconc=bootstrp(10000,@mean,cfpconcinbin); %sample size 10,000. sample the 'mean'. use 'cfpconcinbin' as input distribution
        currentbootmeanyfpconc=bootstrp(10000,@mean,yfpconcinbin);
        binnedbootmeancfpconc=[binnedbootmeancfpconc, std(currentbootmeancfpconc)];  % standard deviation of the sampled mean
        binnedbootmeanyfpconc=[binnedbootmeanyfpconc, std(currentbootmeanyfpconc)];    
        clear cfpconcinbin yfpconcinbin currentbootmeancfpconc currentbootmeanyfpconc
    end
end

figure(12)
clf
plot(phconc(idx),cfpconcuse(idx),'.', 'Color',[0 0.7 1])
hold on
plot(meanbinconc,binnedcfpconc,'.-k','MarkerSize',15)
    
figure(13)
clf
plot(phconc(idx),yfpconcuse(idx),'.','Color',[1 0.7 0]) %[0.9 0.8 0])
%%blubb
%plot(phconc(idx),yfpconcuse(idx),'.','Color',[0.7 0.7 0.7]) %[0.9 0.8 0]) %blubb
hold on
plot(meanbinconc,binnedyfpconc,'.-k','MarkerSize',15)  %blubb
%plot(meanbinconc,binnedyfpconc,'-r','LineWidth',4)  %blubb
if ~isempty(highlightschnitzes)
    for j=1:length(highlightschnitzes)
       ss=highlightschnitzes(j);
       idxschn=find(schnsub==ss);
       phschn=phsub(idxschn);
       yfpconcschn=yfpconcsub(idxschn);
       mycolor=rand(3,1); %blubb
       %if ss==210 | ss==182 | ss==423 %blubb
       %     plot(phschn,yfpconcschn,'.-','MarkerSize',15,'Color',[0 0 0],'LineWidth',3,'MarkerSize',20) % blubb    
       %else %blubb
       %     plot(phschn,yfpconcschn,'.-','MarkerSize',15,'Color',mycolor,'LineWidth',2) % blubb  %blubb
       %end %blubb
       %plot(phschn,yfpconcschn,'.-','MarkerSize',15,'Color',[0.5 0.4 0]) %blubb
    end
end


figure(14)
clf
errorbar(meanbinconc,binnedcfpconc,binnedbootmeancfpconc,'.-','MarkerSize',15,'Color',[0 0.7 1]) %don't divide by mean, at least not in 'noise'case
hold on
errorbar(meanbinconc,binnedyfpconc,binnedbootmeancfpconc,'.-','MarkerSize',15,'Color' ,[1 0.7 0])%[0.9 0.8 0])
grid on

% ************************************************************************************
% ************************************************************************************

%% CELLULAR LENGTH

% Total length increase
phlength=[];
totlength=[];
timedata=[];
schn=[];
frames=[];
frames_unique=[];
%highlightschnitzes=[111 200 444];
%highlightcolors=[0.6 0 0; 0.6 0 0; 0.6 0 0]; % must have the same #columns as length of highlightschnitzes
RESTRICTTOFLUOIMAGES=0;  % only take data which appears at fluorphasefield
fluoframesall='Y_frames_all'; % concentration/fluo image time points. late cells may have one fluorate datapoint less
                           % when no fluo image exists, value is 'nan'

myschnitzcells=s_rm_fitTime;

for i=1:length(myschnitzcells)
    s=myschnitzcells(i);
    if s.useForPlot==1
        if s.completeCycle==1 & ~isempty(s.time) & length(s.phase)==length(s.length_fitNew) & length(s.phase)==length(s.time) 
            if ~RESTRICTTOFLUOIMAGES %take all data
                phlength=[phlength; s.phase'];
                totlength=[totlength; s.length_fitNew']; %smoothing?            
                %totlength=[totlength; s.length_fitNew'/mean(s.length_fitNew)];
                             % normalize by average
                %totlength=[totlength; s.length_fitNew'/(s.length_fitNew(1))];
                            % normalize by initial value
                reltime=s.time-s.time(1);
                timedata=[timedata; reltime'];
                schnrep=zeros(1,length(s.phase))+i;
                schn=[schn; schnrep'];
                framesrep=zeros(1,length(s.phase))+length(s.frames);
                frames=[frames; framesrep'];
                frames_unique=[frames_unique; length(s.frames)];
            else %take only data at fluotimefields
                fluoframeidx=find(~isnan(s.(fluoframesall)));
                phlength=[phlength; s.phase(fluoframeidx)'];
                totlength=[totlength; s.length_fitNew(fluoframeidx)']; %smoothing?            
                %totlength=[totlength; s.length_fitNew(fluoframeidx)'/mean(s.length_fitNew(fluoframeidx))];
                            % normalize by average
                %totlength=[totlength; s.length_fitNew(fluoframeidx)'/(s.length_fitNew(fluoframeidx(1)))];
                            % normalize by initial value
                reltime=s.time(fluoframeidx)-s.time(1);
                timedata=[timedata; reltime'];
                schnrep=zeros(1,length(s.phase(fluoframeidx)))+i;
                schn=[schn; schnrep'];
                framesrep=zeros(1,length(s.phase(fluoframeidx)))+length(s.frames);
                frames=[frames; framesrep'];
                frames_unique=[frames_unique; length(s.frames)];
            end   
        end
    end
end

% normalize and save unnormalized vector
totlengthabs=totlength;
totlength=totlength/mean(totlength);

% binning LENGTH
%idx=find(frames<61 & frames>54);
idx=find(frames<665 & frames>0); %loose conditions like frames>0 & frames<200 include all cells
%plot(phlength(idx),totlength(idx),'.r')
phlengthsub=phlength(idx);
totlengthsub=totlength(idx);


% do binning
numbins=15;
binborders=[0:1/numbins:0.98];
binborders=[binborders,1];
binnedtotlength=[];
meanbinlength=[0.5/numbins:1/numbins:0.98];
binnedstdtotlength=[];
% std error: don't use. see above.
binnederrortotlength=[];
% bootstrapping of mean. use. see above.
binnedbootmeanLength=[];

for i=1:numbins
    idx2=find(phlengthsub>=binborders(i) & phlengthsub<binborders(i+1));
    lengthinbin=totlengthsub(idx2);
    binnedtotlength=[binnedtotlength,mean(lengthinbin)];
    binnedstdtotlength=[binnedstdtotlength,std(lengthinbin)];
    % std erorr: don't use. see above.
    binnederrortotlength=[binnederrortotlength,std(lengthinbin)/sqrt(length(lengthinbin))];    
    % bootstrapping. use. see above. (description for prod. rates)
    currentbootmeanLength=bootstrp(10000,@mean,lengthinbin); %sample size 10,000. sample the 'mean'. use 'lengthinbin' as input distribution
    binnedbootmeanLength=[binnedbootmeanLength, std(currentbootmeanLength)];  % standard deviation of the sampled mean
    clear lengthinbin currentbootmeanLength
end

% fit exponential to bins
[expo, L0]=DJK_ExponentialFit(meanbinlength,binnedtotlength);
phasecont=0:0.01:1;
lengthcont=L0*2.^(expo*phasecont);

%binned length increase
figure(2)
clf
errorbar(meanbinlength,binnedtotlength,binnedbootmeanLength,'.-','MarkerSize',15,'Color','r')
hold on
plot(phasecont, lengthcont,'k')
grid on


%highlightschnitzes=[444 555 666];
figure(1)
clf
hold on
% ***** plot dots
plot(phlength(idx),totlength(idx),'.','Color','r')
%plot lines
%for j=1:length(unique(schn))
%    schnunique=unique(schn);
%    ss=schnunique(j);
%    idxschn=find(schn==ss);
%    phlengthschn=phlength(idxschn);
%    totlengthschn=totlength(idxschn);
%    plot(phlengthschn,totlengthschn,'-','Color','r')
%end
% *****
plot(meanbinlength,binnedtotlength,'.-k','MarkerSize',15)
if ~isempty(highlightschnitzes)
   % if size(highlightcolors,1)==length(highlightschnitzes);
        for j=1:length(highlightschnitzes)
            ss=highlightschnitzes(j);
            idxschn=find(schn==ss);
            phlengthschn=phlength(idxschn);
            totlengthschn=totlength(idxschn);
            plot(phlengthschn,totlengthschn,'.-','Color',[0.5 0 0])
        end
 %   else
 %       disp('wrong # of highlighted lineage colours')
  %  end
end
    




%% RATES RFP GFP

% ************************************************************************************************
% ** THE SCRIPT FOR RFP (MCHERRY) AND GFP IS OUTDATED AND NOT UPDATED ANY  
% ** MORE. TO ANALYZE MOVIES WITH -r- AND -g- , COPY THE RELEVANT FIELDS IN
% ** THE SCHNITZCELLS FILE AND RENAME THEM TO -y- AND -c- . E.g. dR5 -> dC5
% ** AND G5_mean -> Y5_mean. KEEP TRACK OF WHICH COLOR IS TRANSFERED TO WHICH COLOR.
% ************************************************************************************************

% ph=[];
% rfp=[];
% gfp=[];
% time=[];
% schn=[];
% frames=[];
% frames_unique=[];
% abstime=[];
% 
% myschnitzcells=s_rm_fitTime;
% %myschnitzcells=s_rm_fitTime_birthdiv;
% 
% 
% for i=1:length(myschnitzcells)
%     s=myschnitzcells(i);
%     if s.useForPlot==1
%         if length(s.time_atdR)>0 & length(s.phase2_atdR)==length(s.dR5_sum_dt) & length(s.phase2_atdR)==length(s.time_atdR) &  s.completeCycle==1
%             abstime=[abstime; s.time_atdR'];
%             ph=[ph; s.phase2_atdR']; %BLUBB!!! wrong phase/time maybe! !!! fits to dR5_sum_dt
%             rfp=[rfp; s.dR5'];%[rfp; s.dR5_sum_dt'];%[rfp; s.dR5_sum_dt_s'];
%             gfp=[gfp; s.dG5'];%[gfp; s.dG5_sum_dt'];%[rfp; s.dG5_sum_dt_s'];
%             %rfp=[rfp; s.dR5_sum_dt'];
%             %gfp=[gfp; s.dG5_sum_dt'];
%             %rfp=[rfp; s.dR5_sum_dt'/mean(s.dR5_sum_dt)];
%             %gfp=[gfp; s.dG5_sum_dt'/mean(s.dG5_sum_dt)];
%             reltime=s.time_atdR-s.time(1);
%             time=[time; reltime'];
%             schnrep=zeros(1,length(s.phase2_atdR))+i;
%             schn=[schn; schnrep'];
%             framesrep=zeros(1,length(s.phase2_atdR))+length(s.frames);
%             frames=[frames; framesrep'];
%             frames_unique=[frames_unique; length(s.frames)];
%         end
%     end
% end
% 
% % normalize and save unnormalized vector
% rfpabs=rfp;
% gfpabs=gfp;
% rfp=rfp/mean(rfp);
% gfp=gfp/mean(gfp);
% 
% % get noise values (population mean at each time point subtracted)
% % all cells weighed equally
% timeuniquerate=unique(abstime);
% noiserfp=zeros(size(rfp));
% noisegfp=zeros(size(gfp));
% timeTolerance=0.5; % sometimes images taken a little later/earlier
% meangfprate_vec=[];
% meanrfprate_vec=[];
% for i=1:length(timeuniquerate)
%     time=timeuniquerate(i);
%     idxt=find(abstime>time-timeTolerance & abstime<time+timeTolerance);
%     meanrfp_at_t=mean(rfp(idxt));
%     meangfp_at_t=mean(gfp(idxt));
%     noiserfp(idxt)=rfp(idxt)-meanrfp_at_t;
%     noisegfp(idxt)=gfp(idxt)-meangfp_at_t;
%     meanrfprate_vec=[meanrfprate_vec, meanrfp_at_t];
%     meangfprate_vec=[meangfprate_vec, meangfp_at_t];
%     
% end
% 
% figure(8)
% clf
% plot(timeuniquerate,meangfprate_vec,'.-','MarkerSize',10,'Color',[0 0.7 0]) %don't divide by mean, at least not in 'noise'case
% hold on
% plot(timeuniquerate,meanrfprate_vec,'.-r','MarkerSize',10)
% title('time dependence av rates')
% grid on
% 
% 
% % choose if abs value (normalized) or noise (mean subtracted) is used
% rfpuse=rfp;
% gfpuse=gfp;
% %rfpuse=noiserfp;
% %gfpuse=noisegfp;
% 
% 
% % binning RATES RFP GFP
% framemin=0;
% framemax=890;
% %idx=find(frames<55 & frames>35);
% idx=find(frames<framemax & frames>framemin);
% 
% phsub=ph(idx);
% rfpsub=rfpuse(idx);
% gfpsub=gfpuse(idx);
% 
% 
% % do binning
% numbinsrate=8;
% numbins=numbinsrate;
% binborders=[0:1/numbins:0.98];
% binborders=[binborders,1];
% binnedgfp=[];
% binnedrfp=[];
% meanbin=[0.5/numbins:1/numbins:0.98];
% 
% for i=1:numbins
%     idx2=find(phsub>=binborders(i) & phsub<binborders(i+1));
%     binnedrfp=[binnedrfp,mean(rfpsub(idx2))];
%     binnedgfp=[binnedgfp,mean(gfpsub(idx2))];
% end
% 
% figure(9)
% clf
% plot(ph(idx),rfp(idx),'.r')
% hold on
% plot(meanbin,binnedrfp,'.-k','MarkerSize',15)
% %plot(meanbin,4*binnedrfp,'.-','MarkerSize',10,'Color',[0.5 0.5 0.5])
% %plot(meanbin,6*binnedrfp,'.-','MarkerSize',10,'Color',[0.5 0.5 0.5])
%     
% figure(10)
% clf
% plot(ph(idx),gfp(idx),'.','Color',[0 0.7 0])
% hold on
% plot(meanbin,binnedgfp,'.-k','MarkerSize',15)
% %plot(meanbin,4*binnedgfp,'.-','MarkerSize',10,'Color',[0.5 0.5 0.5])
% %plot(meanbin,6*binnedgfp,'.-','MarkerSize',10,'Color',[0.5 0.5 0.5])
% 
% figure(11)
% clf
% plot(meanbin,binnedrfp,'.-r','MarkerSize',15)
% hold on
% plot(meanbin,binnedgfp,'.-','MarkerSize',15,'Color', [0 0.7 0])
% grid on
% 
% 
% % CONC RFP GFP
% 
% phconc=[];
% rfpconc=[];
% gfpconc=[];
% timeconc=[];
% schnconc=[];
% framesconc=[];
% frames_uniqueconc=[];
% abstimeconc=[];
% 
% 
% 
% 
% 
% for i=1:length(myschnitzcells)
%     s=myschnitzcells(i);
%     if s.useForPlot==1
%         if length(s.R_time)>0 & length(s.phase2_atR)==length(s.R5_mean) & length(s.phase2_atR)==length(s.R_time) &  s.completeCycle==1
%             phconc=[phconc; s.phase2_atR'];
%             %rfpconc=[rfpconc; s.R5_mean'/mean(s.R5_mean)];
%             %gfpconc=[gfpconc; s.G5_mean'/mean(s.G5_mean)];
%             abstimeconc=[abstimeconc; s.R_time'];
%             rfpconc=[rfpconc; s.R5_mean'];
%             gfpconc=[gfpconc; s.G5_mean'];
%             
%             reltimeconc=s.R_time-s.time(1);
%             timeconc=[timeconc; reltime'];
%             schnrepconc=zeros(1,length(s.phase2_atR))+i;
%             schnconc=[schnconc; schnrepconc'];
%             framesrepconc=zeros(1,length(s.phase2_atR))+length(s.frames);
%             framesconc=[framesconc; framesrepconc'];
%             frames_uniqueconc=[frames_uniqueconc; length(s.frames)];
%         end
%     end
% end
% 
% % normalize and save unnormalized vector
% rfpconcabs=rfpconc;
% gfpconcabs=gfpconc;
% rfpconc=rfpconc/mean(rfpconc);
% gfpconc=gfpconc/mean(gfpconc);
% 
% % get noise values (population mean at each time point subtracted)
% % all cells weighed equally
% timeunique=unique(abstimeconc);
% noiserfpconc=zeros(size(rfpconc));
% noisegfpconc=zeros(size(gfpconc));
% timeTolerance=0.5; % sometimes images taken a little later/earlier
% meangfp_vec=[];
% meanrfp_vec=[];
% for i=1:length(timeunique)
%     time=timeunique(i);
%     idxt=find(abstimeconc>time-timeTolerance & abstimeconc<time+timeTolerance);
%     meanrfpconc_at_t=mean(rfpconc(idxt));
%     meangfpconc_at_t=mean(gfpconc(idxt));
%     noiserfpconc(idxt)=rfpconc(idxt)-meanrfpconc_at_t;
%     noisegfpconc(idxt)=gfpconc(idxt)-meangfpconc_at_t;
%     meanrfp_vec=[meanrfp_vec, meanrfpconc_at_t];
%     meangfp_vec=[meangfp_vec, meangfpconc_at_t];
%     
% end
% 
% figure(15)
% clf
% plot(timeunique,meangfp_vec,'.-','MarkerSize',10,'Color',[0 0.7 0]) %don't divide by mean, at least not in 'noise'case
% hold on
% plot(timeunique,meanrfp_vec,'.-r','MarkerSize',10)
% title('time dependence av conc R5 etc')
% grid on
% 
% 
% % binning conc RFP GFP
% 
% % choose if abs value (normalized) or noise (mean subtracted) is used
% %rfpconcuse=rfpconc;
% %gfpconcuse=gfpconc;
% rfpconcuse=noiserfpconc;
% gfpconcuse=noisegfpconc;
% 
% %idx=find(frames<49 & frames>0); % choose framerange
% %idx=find(frames<59 & frames>45) % CHOSEN ABOVE! BLUBB
% idx=find(frames<framemax & frames>framemin); % CHOSEN ABOVE!
% 
% phconcsub=phconc(idx);
% gfpconcsub=gfpconcuse(idx);
% rfpconcsub=rfpconcuse(idx);
% 
% % choose if binning is supposed to be periodic (conc is continuous!)
% periodic=0; % =1; =0;
% 
% if periodic
%     % do binning periodic (pt at =0 and =1 are identical)
%     numbins=7; %actually only (n-1) independent bins
%     binborders=[0.5/(numbins-1):1/(numbins-1):0.98];
%     meanbin=[0:1/(numbins-1):0.98];
%     meanbin=[meanbin, 1];
%     binnedrfpconc=[];
%     binnedgfpconc=[];
% 
%     for i=1:(numbins-1)
%         if i==1 % first=last idx, since periodic binning
%             idx2=find(phconcsub<binborders(i) | phconcsub>=binborders(end));
%             binnedgfpconc=[binnedgfpconc,mean(gfpconcsub(idx2))];
%             binnedrfpconc=[binnedrfpconc,mean(rfpconcsub(idx2))];
%         else
%             idx2=find(phconcsub>=binborders(i-1) & phconcsub<binborders(i));
%             binnedgfpconc=[binnedgfpconc,mean(gfpconcsub(idx2))];
%             binnedrfpconc=[binnedrfpconc,mean(rfpconcsub(idx2))];
%         end
%     end
%     % at last binned data (periodic)
%      binnedgfpconc=[binnedgfpconc,binnedgfpconc(1)];
%      binnedrfpconc=[binnedrfpconc,binnedrfpconc(1)];
% else
%     % do binning not=periodic
%     numbins=numbinsrate;
%     binborders=[0:1/numbins:0.98];
%     binborders=[binborders,1];
%     binnedrfpconc=[];
%     binnedgfpconc=[];
%     meanbin=[0.5/numbins:1/numbins:0.98];
% 
%     for i=1:numbins
%         idx2=find(phconcsub>=binborders(i) & phconcsub<binborders(i+1));
%         binnedgfpconc=[binnedgfpconc,mean(gfpconcsub(idx2))];
%         binnedrfpconc=[binnedrfpconc,mean(rfpconcsub(idx2))];
%     end
% end
% 
% figure(12)
% clf
% plot(phconc(idx),gfpconcuse(idx),'.', 'Color',[0 0.7 0])
% hold on
% plot(meanbin,binnedgfpconc,'.-k','MarkerSize',15)
%     
% figure(13)
% clf
% plot(phconc(idx),rfpconcuse(idx),'.r')
% hold on
% plot(meanbin,binnedrfpconc,'.-k','MarkerSize',15)
% 
% figure(14)
% clf
% plot(meanbin,binnedgfpconc,'.-','MarkerSize',15,'Color',[0 0.7 0]) %don't divide by mean, at least not in 'noise'case
% hold on
% plot(meanbin,binnedrfpconc,'.-r','MarkerSize',15)
% grid on
% 


