% This script is used to determine the bleaching rates of GFP and mCherry
% proteins.
% --------------------------------------------------------------------
% It is written specificially for analysis of the experiment from
% 2014-03-27 (i.e. naming is specific to the labels) but could easily be
% altered for other experiments.
% Assumed datasets: (a) rpoS-GFP, L31-mCherry, (b) icd-GFP, L31-mCherry.
%
% The experimental data: Bacteria were treated with a high conc of
% antibiotics so that protein production stopped. Subsequently they were 
% fluo-imaged many times (70x in 2014-03-27) with a standard time-lapse
% journal. Fluo concentration (GFP and mCherry) of cells follows an
% exponential decrease. Data was analyzed with schnitzcells software.
%
% In this script the schnitzcells data is used to determine the exponential
% decay constants (bleaching per illumination) and bleaching per second (by
% using the illlumination times). It is checked that different
% Fluo-quantities (G6mean,G5mean,G5sum) does not alter anything. [Note
% G5sum is sometimes names G5total].
% Using rpos-gfp and icd-gfp (very different signal strength) allows to
% determine the influence of initial fluorophore concentration.
%
% ********* STEPS *********
% ** run each cell seperately **
% 1) load schnitzcells data and specify illumination times
% 2) extract fluorescence-over-time matrices from schnitzcell vectors
% 3) (optional) test whether low-illuminated control cells had constant
%     fluroescence
% 4) determine bleaching rate per illumination event. convert to
%     bleaching rate per second. plot bleaching rates and curves. 
%    check dependence on fluorescence parameter and GFP label
%    [add on: determine noise around fitted bleaching curve as a proxy for
%    measurement noise (assuming that bleaching is a quasi-deterministic
%    decreasing exponential)]
% 5) combine data and final bleach value


%% -----------------------------------------------------------------------
% (1) LOAD SCHNITZCELLS. SPECIFY PATH AND ILLUMINATION TIME
% ------------------------------------------------------------------------
% path
myrootdir='\\biofysicasrv\Users2\Walker\Experiments\Technical\2014-03-27\';
allfiles=dir([myrootdir 'schnitz*.mat']);
for i=1:length(allfiles)
    load([myrootdir allfiles(i).name] )
end

% illumination times [ms]
illumTimeGFP_icd=300;
illumTimeGFP_rpos=300;
illumTimemCherry_L31=150;


clear allfiles i
%% --------------------------------------------------------------------------
% (2) EXTRACT FLUORESCENCE-OVER-TIME MATRICES FROM SCHNITZCELLS FILES
% Specific to 2014-03-27 schnitzcells names
% ---------------------------------------------------------------------------

% structure, e.g. G6_mean_icdL31
%            schnitz1  -  schnitz2  -  schnitz3  -  schnitz4  -  ...
% 1st row:     G6_t0        G6_t0        G6_t0        G6_t0
% 2nd row:     G6_t1        G6_t1        G6_t1        G6_t1
% 3rd row:      ...         ....          ...          ...

% ------------------------------
% ICD & L31.
% ------------------------------
% *** define dataset ***
myschnitzcells=schnitz_icd_L31;
mysize=[length(myschnitzcells(1).G6_mean),length(myschnitzcells)];

% *** allocate empty matrices ***
% GFP
G6mean_icdL31=zeros(mysize);
G5mean_icdL31=zeros(mysize);
G5total_icdL31=zeros(mysize);
% mCherry
R6mean_icdL31=zeros(mysize);
R5mean_icdL31=zeros(mysize);
R5total_icdL31=zeros(mysize);

% *** fill matrices with data ***
for i=1:length(myschnitzcells)
    G6mean_icdL31(:,i)=myschnitzcells(i).G6_mean';
    G5mean_icdL31(:,i)=myschnitzcells(i).G5_mean';
    G5total_icdL31(:,i)=myschnitzcells(i).G5_sum';
    
    R6mean_icdL31(:,i)=myschnitzcells(i).R6_mean';
    R5mean_icdL31(:,i)=myschnitzcells(i).R5_mean';
    R5total_icdL31(:,i)=myschnitzcells(i).R5_sum';
end


% ------------------------------
% ICD & L31. LOW ILLUMINATION CONTROL: Only G/R5total examined
% ------------------------------
% *** define dataset ***
myschnitzcells=schnitz_icd_L31_lowfrequ;
mysize=[length(myschnitzcells(1).G6_mean),length(myschnitzcells)];

% *** allocate empty matrices ***
% GFP
G5total_icdL31_lowfrequ=zeros(mysize);
% mCherry
R5total_icdL31_lowfrequ=zeros(mysize);

% *** fill matrices with data ***
for i=1:length(myschnitzcells)
    G5total_icdL31_lowfrequ(:,i)=myschnitzcells(i).G5_sum';
    R5total_icdL31_lowfrequ(:,i)=myschnitzcells(i).R5_sum';
end



% ------------------------------
% RPOS & L31.
% ------------------------------
% *** define dataset ***
myschnitzcells=schnitz_rpos_L31;
mysize=[length(myschnitzcells(1).G6_mean),length(myschnitzcells)];

% *** allocate empty matrices ***
% GFP
G6mean_rposL31=zeros(mysize);
G5mean_rposL31=zeros(mysize);
G5total_rposL31=zeros(mysize);
% mCherry
R6mean_rposL31=zeros(mysize);
R5mean_rposL31=zeros(mysize);
R5total_rposL31=zeros(mysize);

% *** fill matrices with data ***
for i=1:length(myschnitzcells)
    G6mean_rposL31(:,i)=myschnitzcells(i).G6_mean';
    G5mean_rposL31(:,i)=myschnitzcells(i).G5_mean';
    G5total_rposL31(:,i)=myschnitzcells(i).G5_sum';
    
    R6mean_rposL31(:,i)=myschnitzcells(i).R6_mean';
    R5mean_rposL31(:,i)=myschnitzcells(i).R5_mean';
    R5total_rposL31(:,i)=myschnitzcells(i).R5_sum';
end


% ------------------------------
% RPOS & L31. LOW ILLUMINATION CONTROL: Only G/R5total examined
% ------------------------------
% *** define dataset ***
myschnitzcells=schnitz_rpos_L31_lowfrequ;
mysize=[length(myschnitzcells(1).G6_mean),length(myschnitzcells)];

% *** allocate empty matrices ***
% GFP
G5total_rposL31_lowfrequ=zeros(mysize);
% mCherry
R5total_rposL31_lowfrequ=zeros(mysize);

% *** fill matrices with data ***
for i=1:length(myschnitzcells)
    G5total_rposL31_lowfrequ(:,i)=myschnitzcells(i).G5_sum';
    R5total_rposL31_lowfrequ(:,i)=myschnitzcells(i).R5_sum';
end



% some illustrative figures of the traces
figure(1)
clf
hold on
plot(G5total_icdL31,'LineWidth',2)
xlabel('x''th illumination step')
ylabel('total GFP. icd')
figure(2)
clf
hold on
plot(R5total_icdL31,'LineWidth',2)
xlabel('x''th illumination step')
ylabel('total mCherry. l31 (icd strain)')
figure(3)
clf
hold on
plot(G5total_rposL31,'LineWidth',2)
xlabel('x''th illumination step')
ylabel('total GFP. rpos')
figure(4)
clf
hold on
plot(R5total_rposL31,'LineWidth',2)
xlabel('x''th illumination step')
ylabel('total mCherry. l31 (rpos strain)')

%% -----------------------------------------------------------------------
% (3) CHECK CONSTANT FLUORESCENCE OF LOW ILLUMINATION CONTROL BACTERIA
% ------------------------------------------------------------------------
% optional

% ------------------------------
% ICD & L31.
% ------------------------------
figure(5)
clf
hold on
subplot(1,2,1)
plot(G5total_icdL31_lowfrequ/mean(mean(G5total_icdL31_lowfrequ)),'LineWidth',2)
xlabel('x''th illumination')
ylabel('total GFP (G5sum) normalized by mean')
title('ICD-GFP.  20x less frequently illuminated. single cells')
subplot(1,2,2)
plot(R5total_icdL31_lowfrequ/mean(mean(R5total_icdL31_lowfrequ)),'LineWidth',2)
xlabel('x''th illumination')
ylabel('total mCherry (R5sum) normalized by mean')
title('L31-mCherry.  20x less frequently illuminated. single cells')

% ------------------------------
% RPOS & L31.
% ------------------------------
figure(6)
clf
hold on
subplot(1,2,1)
plot(G5total_rposL31_lowfrequ/mean(mean(G5total_rposL31_lowfrequ)),'LineWidth',2)
xlabel('x''th illumination')
ylabel('total GFP (G5sum) normalized by mean')
title('RPOS-GFP.  20x less frequently illuminated. single cells')
subplot(1,2,2)
plot(R5total_rposL31_lowfrequ/mean(mean(R5total_rposL31_lowfrequ)),'LineWidth',2)
xlabel('x''th illumination')
ylabel('total mCherry (R5sum) normalized by mean')
title('L31-mCherry.  20x less frequently illuminated. single cells')




%% ----------------------------------------------------------------------------
% (4) DETERMINE BLEACH RATE PER ILLUMINATION EVENT. CONVERT TO BLEACH RATE
% PER SECOND. PLOT.
% RUN ONCE PER DATASET AND CHANGE THE CHOICE OF YOUR DATASET
% ----------------------------------------------------------------------------
% Analysis is performed for all possible fluroescence measures: G6mean,G5mean,G5total. Same for
% R5/6.
% Only the bleach rates determined with total fluroescence (G5sum,R5sum
% resp. "total") is stored for further analysis because it will show that
% choice of fluo-variable is irrelevant.

% This script is largely copied from 'CompleteScriptMaturationTimes'
% (subsection Bleaching)

% plot fitted bleach curve for each cell.
PLOTALLTRACES=1;

% ********** CHOOSE DATASET *****************
MYDATASET='ICD_L31';
% MYDATASET='RPOS_L31';

% *************************************************************************
if (strcmp(MYDATASET,'ICD_L31')==1)
    G6mean=G6mean_icdL31;
    G5mean=G5mean_icdL31;
    G5total=G5total_icdL31;
    R6mean=R6mean_icdL31;
    R5mean=R5mean_icdL31;
    R5total=R5total_icdL31;
    mytitle1='ICD_L31: ';
    illumTimeGFP=illumTimeGFP_icd;
    illumTimemCherry=illumTimemCherry_L31;
    % output variable names. They will refer to results of G5total/R5total
    bleach_per_sec_GFP='bleach_per_sec_GFP_icd'; % exponential bleach rate per second
    bleach_per_sec_mCherry='bleach_per_sec_mCherry_L31_1'; % cannot use same name for both strains -> _1,_2
    bleach_per_illum_GFP='bleach_per_illum_GFP_icd'; % bleach factor for one illumination with this settings (illum time)
    bleach_per_illum_mCherry='bleach_per_illum_mCherry_L31_1';
    % ************************
elseif (strcmp(MYDATASET,'RPOS_L31')==1)
    G6mean=G6mean_rposL31;
    G5mean=G5mean_rposL31;
    G5total=G5total_rposL31;
    R6mean=R6mean_rposL31;
    R5mean=R5mean_rposL31;
    R5total=R5total_rposL31;
    mytitle1='RPOS_L31: ';
    illumTimeGFP=illumTimeGFP_rpos;
    illumTimemCherry=illumTimemCherry_L31;
    % output variable names. They will refer to results of G5total/R5total
    bleach_per_sec_GFP='bleach_per_sec_GFP_rpos'; % exponential bleach rate per second
    bleach_per_sec_mCherry='bleach_per_sec_mCherry_L31_2'; % cannot use same name for both strains -> _1,_2
    bleach_per_illum_GFP='bleach_per_illum_GFP_rpos'; % bleach factor for one illumination with this settings (illum time)
    bleach_per_illum_mCherry='bleach_per_illum_mCherry_L31_2';
else
    error('Unknown dataset chosen.')
end

% *************************************************************************
% from here on the determination of bleaching rates is for general vector
% names.
% *************************************************************************

% *** create dummy time vector ***
% only works if same # datapoints for each vector
xvec=[1:size(G6mean,1)]';

% *** define a long list of colors for plotting traces ***
myColor=[1 1 0.9; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];

% *** create vectors for bleach rates ***
G6meanbleachrates=zeros(size(G6mean,2),1);
G5meanbleachrates=zeros(size(G6mean,2),1);
G5totalbleachrates=zeros(size(G6mean,2),1);
R6meanbleachrates=zeros(size(G6mean,2),1);
R5meanbleachrates=zeros(size(G6mean,2),1);
R5totalbleachrates=zeros(size(G6mean,2),1);

% [add on: *** create vectors for noise in concentration ***]
G6meannoise=zeros(size(G6mean,2),1);
G5meannoise=zeros(size(G6mean,2),1);
G5totalnoise=zeros(size(G6mean,2),1);
R6meannoise=zeros(size(G6mean,2),1);
R5meannoise=zeros(size(G6mean,2),1);
R5totalnoise=zeros(size(G6mean,2),1);


% *************************************************************************
% *** G6mean ***
% loop over all schnitzes
for i=1:size(G6mean,2)
    % *** adjust ***
    fluovec=G6mean(:,i);
    mytitle2='G6mean. ';
    % ***************
    
    % fit an exponential bleaching function (adapted from curve fitting
    % session & maturation time script)
    ok2 = isfinite(xvec) & isfinite(fluovec);
    if ~all( ok2 )
        warning( 'Bleach factor calculation: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
    end
    % suitable start paramters for [a b]
    st2= [ 400000 -0.01];
    ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
    % Fit this model using new data
    cfbleach = fit(xvec(ok2),fluovec(ok2),ft2,'Startpoint',st2);

    % get fitted line
    a=cfbleach.a; bleachrate=-cfbleach.b;
    timeaxisbleach=[1:max(xvec)+2];
    yaxisfitbleach=a*exp(-bleachrate*timeaxisbleach);
    if PLOTALLTRACES
        %figure(3)
        %if i==1
        %    clf
        %    hold on
        %    set(gcf,'WindowStyle','docked')
        %end
        %h1=plot(xvec,fluovec,'.','MarkerSize',15,'Color',myColor(i,:));
        %h2=plot(timeaxisbleach,yaxisfitbleach,'LineWidth',2,'Color',myColor(i,:));
        figure
        clf 
        hold on
        set(gcf,'WindowStyle','docked')
        h1=plot(xvec,fluovec,'.','MarkerSize',15);
        h2=plot(timeaxisbleach,yaxisfitbleach,'r','LineWidth',2);
        xlabel('pseudo time (illum steps)')
        ylabel('fluorescence')
        legend(h2,['b=-' num2str(bleachrate)])
        title(['Bleach curve: ' mytitle1 mytitle2 ' schnitz nr ' num2str(i)],'Interpreter','None');
    end
    
    % store all bleachrates into one vector
    G6meanbleachrates(i)=bleachrate;
    
    %calculate deviations from fit (presumably mostly measurement error /
    %lamp fluctuations
    residuals=fluovec-a*exp(-bleachrate*xvec);
    % get noise: stddev, then normalize by average fluo-value
    noisefluo=std(residuals)/mean(fluovec);
    G6meannoise(i)=noisefluo;
    
end

% *************************************************************************
% *** G5mean ***
% loop over all schnitzes
for i=1:size(G5mean,2)
    % *** adjust ***
    fluovec=G5mean(:,i);
    mytitle2='G5mean. ';
    % ***************
    
    % fit an exponential bleaching function (adapted from curve fitting
    % session & maturation time script)
    ok2 = isfinite(xvec) & isfinite(fluovec);
    if ~all( ok2 )
        warning( 'Bleach factor calculation: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
    end
    % suitable start paramters for [a b]
    st2= [ 400000 -0.01];
    ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
    % Fit this model using new data
    cfbleach = fit(xvec(ok2),fluovec(ok2),ft2,'Startpoint',st2);

    % get fitted line
    a=cfbleach.a; bleachrate=-cfbleach.b;
    timeaxisbleach=[1:max(xvec)+2];
    yaxisfitbleach=a*exp(-bleachrate*timeaxisbleach);
    if PLOTALLTRACES
        figure
        clf
        hold on
        set(gcf,'WindowStyle','docked')
        h1=plot(xvec,fluovec,'.','MarkerSize',15);
        h2=plot(timeaxisbleach,yaxisfitbleach,'r','LineWidth',2);
        xlabel('pseudo time (illum steps)')
        ylabel('fluorescence')
        legend(h2,['b=-' num2str(bleachrate)])
        title(['Bleach curve: ' mytitle1 mytitle2 ' schnitz nr ' num2str(i)],'Interpreter','None');
    end
    
    % store all bleachrates into one vector
    G5meanbleachrates(i)=bleachrate;
    
    %calculate deviations from fit (presumably mostly measurement error /
    %lamp fluctuations
    residuals=fluovec-a*exp(-bleachrate*xvec);
    % get noise: stddev, then normalize by average fluo-value
    noisefluo=std(residuals)/mean(fluovec);
    G5meannoise(i)=noisefluo;
    
end


% *************************************************************************
% *** G5total ***
% loop over all schnitzes
for i=1:size(G5total,2)
    % *** adjust ***
    fluovec=G5total(:,i);
    mytitle2='G5total. ';
    % ***************
    
    % fit an exponential bleaching function (adapted from curve fitting
    % session & maturation time script)
    ok2 = isfinite(xvec) & isfinite(fluovec);
    if ~all( ok2 )
        warning( 'Bleach factor calculation: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
    end
    % suitable start paramters for [a b]
    st2= [ 400000 -0.01];
    ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
    % Fit this model using new data
    cfbleach = fit(xvec(ok2),fluovec(ok2),ft2,'Startpoint',st2);

    % get fitted line
    a=cfbleach.a; bleachrate=-cfbleach.b;
    timeaxisbleach=[1:max(xvec)+2];
    yaxisfitbleach=a*exp(-bleachrate*timeaxisbleach);
    if PLOTALLTRACES
        figure
        clf
        hold on
        set(gcf,'WindowStyle','docked')
        h1=plot(xvec,fluovec,'.','MarkerSize',15);
        h2=plot(timeaxisbleach,yaxisfitbleach,'r','LineWidth',2);
        xlabel('pseudo time (illum steps)')
        ylabel('fluorescence')
        legend(h2,['b=-' num2str(bleachrate)])
        title(['Bleach curve: ' mytitle1 mytitle2 ' schnitz nr ' num2str(i)],'Interpreter','None');
    end
    
    % store all bleachrates into one vector
    G5totalbleachrates(i)=bleachrate;

    %calculate deviations from fit (presumably mostly measurement error /
    %lamp fluctuations
    residuals=fluovec-a*exp(-bleachrate*xvec);
    % get noise: stddev, then normalize by average fluo-value
    noisefluo=std(residuals)/mean(fluovec);
    G5totalnoise(i)=noisefluo;
    
end



% *************************************************************************
% *** R6mean ***
% loop over all schnitzes
for i=1:size(R6mean,2)
    % *** adjust ***
    fluovec=R6mean(:,i);
    mytitle2='R6mean. ';
    % ***************
    
    % fit an exponential bleaching function (adapted from curve fitting
    % session & maturation time script)
    ok2 = isfinite(xvec) & isfinite(fluovec);
    if ~all( ok2 )
        warning( 'Bleach factor calculation: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
    end
    % suitable start paramters for [a b]
    st2= [ 400000 -0.01];
    ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
    % Fit this model using new data
    cfbleach = fit(xvec(ok2),fluovec(ok2),ft2,'Startpoint',st2);

    % get fitted line
    a=cfbleach.a; bleachrate=-cfbleach.b;
    timeaxisbleach=[1:max(xvec)+2];
    yaxisfitbleach=a*exp(-bleachrate*timeaxisbleach);
    if PLOTALLTRACES
        figure
        clf
        hold on
        set(gcf,'WindowStyle','docked')
        h1=plot(xvec,fluovec,'.','MarkerSize',15);
        h2=plot(timeaxisbleach,yaxisfitbleach,'r','LineWidth',2);
        xlabel('pseudo time (illum steps)')
        ylabel('fluorescence')
        legend(h2,['b=-' num2str(bleachrate)])
        title(['Bleach curve: ' mytitle1 mytitle2 ' schnitz nr ' num2str(i)],'Interpreter','None');
    end
    
    % store all bleachrates into one vector
    R6meanbleachrates(i)=bleachrate;

    %calculate deviations from fit (presumably mostly measurement error /
    %lamp fluctuations
    residuals=fluovec-a*exp(-bleachrate*xvec);
    % get noise: stddev, then normalize by average fluo-value
    noisefluo=std(residuals)/mean(fluovec);
    R6meannoise(i)=noisefluo;
    
end

% *************************************************************************
% *** R5mean ***
% loop over all schnitzes
for i=1:size(R5mean,2)
    % *** adjust ***
    fluovec=R5mean(:,i);
    mytitle2='R5mean. ';
    % ***************
    
    % fit an exponential bleaching function (adapted from curve fitting
    % session & maturation time script)
    ok2 = isfinite(xvec) & isfinite(fluovec);
    if ~all( ok2 )
        warning( 'Bleach factor calculation: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
    end
    % suitable start paramters for [a b]
    st2= [ 400000 -0.01];
    ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
    % Fit this model using new data
    cfbleach = fit(xvec(ok2),fluovec(ok2),ft2,'Startpoint',st2);

    % get fitted line
    a=cfbleach.a; bleachrate=-cfbleach.b;
    timeaxisbleach=[1:max(xvec)+2];
    yaxisfitbleach=a*exp(-bleachrate*timeaxisbleach);
    if PLOTALLTRACES
        figure
        clf
        hold on
        set(gcf,'WindowStyle','docked')
        h1=plot(xvec,fluovec,'.','MarkerSize',15);
        h2=plot(timeaxisbleach,yaxisfitbleach,'r','LineWidth',2);
        xlabel('pseudo time (illum steps)')
        ylabel('fluorescence')
        legend(h2,['b=-' num2str(bleachrate)])
        title(['Bleach curve: ' mytitle1 mytitle2 ' schnitz nr ' num2str(i)],'Interpreter','None');
    end
    
    % store all bleachrates into one vector
    R5meanbleachrates(i)=bleachrate;

    %calculate deviations from fit (presumably mostly measurement error /
    %lamp fluctuations
    residuals=fluovec-a*exp(-bleachrate*xvec);
    % get noise: stddev, then normalize by average fluo-value
    noisefluo=std(residuals)/mean(fluovec);
    R5meannoise(i)=noisefluo;
    
end


% *************************************************************************
% *** R5total ***
% loop over all schnitzes
for i=1:size(R5total,2)
    % *** adjust ***
    fluovec=R5total(:,i);
    mytitle2='R5total. ';
    % ***************
    
    % fit an exponential bleaching function (adapted from curve fitting
    % session & maturation time script)
    ok2 = isfinite(xvec) & isfinite(fluovec);
    if ~all( ok2 )
        warning( 'Bleach factor calculation: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
    end
    % suitable start paramters for [a b]
    st2= [ 400000 -0.01];
    ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
    % Fit this model using new data
    cfbleach = fit(xvec(ok2),fluovec(ok2),ft2,'Startpoint',st2);

    % get fitted line
    a=cfbleach.a; bleachrate=-cfbleach.b;
    timeaxisbleach=[1:max(xvec)+2];
    yaxisfitbleach=a*exp(-bleachrate*timeaxisbleach);
    if PLOTALLTRACES
        figure
        clf
        hold on
        set(gcf,'WindowStyle','docked')
        h1=plot(xvec,fluovec,'.','MarkerSize',15);
        h2=plot(timeaxisbleach,yaxisfitbleach,'r','LineWidth',2);
        xlabel('pseudo time (illum steps)')
        ylabel('fluorescence')
        legend(h2,['b=-' num2str(bleachrate)])
        title(['Bleach curve: ' mytitle1 mytitle2 ' schnitz nr ' num2str(i)],'Interpreter','None');
    end
    
    % store all bleachrates into one vector
    R5totalbleachrates(i)=bleachrate;

    %calculate deviations from fit (presumably mostly measurement error /
    %lamp fluctuations
    residuals=fluovec-a*exp(-bleachrate*xvec);
    % get noise: stddev, then normalize by average fluo-value
    noisefluo=std(residuals)/mean(fluovec);
    R5totalnoise(i)=noisefluo;
    
end


% *************************************************************************
% convert "bleach rate per illum event (b_illum)" (so far calculated) to 
% "bleach rate per second (b_sec)" (general)
% and "bleach factor per illum event b_factor" (specific for this illumination
% time
% *************************************************************************
%
% exp(-b_illum) = exp(-b_sec*time) = b_factor          [time=illum time in sec]
% -> b_sec=b_illum/time
% -> b_factor=exp(-b_illum)

% G5mean
b_sec_G5mean=G5meanbleachrates./illumTimeGFP*1000;  % 1000 for ms->sec conversion
b_factor_G5mean=exp(-G5meanbleachrates);
% G6mean
b_sec_G6mean=G6meanbleachrates./illumTimeGFP*1000;
b_factor_G6mean=exp(-G6meanbleachrates);
% G5total
b_sec_G5total=G5totalbleachrates./illumTimeGFP*1000;
b_factor_G5total=exp(-G5totalbleachrates);
% R5mean
b_sec_R5mean=R5meanbleachrates./illumTimemCherry*1000;
b_factor_R5mean=exp(-R5meanbleachrates);
% R6mean
b_sec_R6mean=R6meanbleachrates./illumTimemCherry*1000;
b_factor_R6mean=exp(-R6meanbleachrates);
% R5total
b_sec_R5total=R5totalbleachrates./illumTimemCherry*1000;
b_factor_R5total=exp(-R5totalbleachrates);


% *************************************************************************
% plot "bleach rate per second" for different fluo-parameters (G5/G6 etc)

avGFPbleach=[mean(b_sec_G6mean) mean(b_sec_G5mean) mean(b_sec_G5total)];
stdGFPbleach=[std(b_sec_G6mean) std(b_sec_G5mean) std(b_sec_G5total)];
avmCherrybleach=[mean(b_sec_R6mean) mean(b_sec_R5mean) mean(b_sec_R5total)];
stdmCherrybleach=[std(b_sec_R6mean) std(b_sec_R5mean) std(b_sec_R5total)];

% GFP
figure(5)
clf
hold on
bar(avGFPbleach, 'FaceColor', [0.3 1 0.3])  % Lighter so error bars show up
errorbar(avGFPbleach,stdGFPbleach, 'ks','LineWidth',2);            % Error bars use black squares
set(gca, 'XTick', 1:3, 'XTickLabel', {'G6mean' 'G5mean'  'G5total'}) % Set ticks and tick labels
ylabel('bleachrate (mean,std)')
title([mytitle1 'GFP. exp bleach rate per sec'],'Interpreter','None')
legend('Mean (SD error bars)', 'Location', 'Northwest') % Put in lower right
box on                                         % Force box around axes
hold off


% mCherry
figure(4)
clf
hold on
bar(avmCherrybleach, 'FaceColor', [1 0.3 0.3])  % Lighter so error bars show up
errorbar(avmCherrybleach,stdmCherrybleach, 'ks','LineWidth',2);            % Error bars use black squares
set(gca, 'XTick', 1:3, 'XTickLabel', {'R6mean' 'R5mean'  'R5total'}) % Set ticks and tick labels
ylabel('bleachrate (mean,std)')
title([mytitle1 'mCherry. exp bleach rate per sec'],'Interpreter','None')
legend('Mean (SD error bars)', 'Location', 'Northwest') % Put in lower right
box on                                         % Force box around axes
hold off


% *************************************************************************
% store bleach rates and factors in a general vector name (which will not
% be overwritten)
% *************************************************************************
% vectors. one value for each schnitz

eval([bleach_per_sec_GFP ' = b_sec_G5total;']) % exponential bleach rate per sec
eval([bleach_per_sec_mCherry ' = b_sec_R5total;'])
eval([bleach_per_illum_GFP ' = b_factor_G5total;']) % bleach factor for given illum time
eval([bleach_per_illum_mCherry ' = b_factor_R5total;'])


%% ----------------------------------------------------------------------------
% (5) MERGE ALL DATA. WRITTEN FOR SPECIFIC VARIABLE NAMES
% ----------------------------------------------------------------------------
% maybe useful to make script more flexible!

%gfp
mean_bleach_per_sec_GFP_rpos=mean(bleach_per_sec_GFP_rpos)
std_bleach_per_sec_GFP_rpos=std(bleach_per_sec_GFP_rpos)
mean_bleach_per_sec_GFP_icd=mean(bleach_per_sec_GFP_icd)
std_bleach_per_sec_GFP_icd=std(bleach_per_sec_GFP_icd)
%merge gfp all together
mean_bleach_per_sec_GFP_merge=mean([ bleach_per_sec_GFP_icd; bleach_per_sec_GFP_rpos])
std_bleach_per_sec_GFP_merge=std([ bleach_per_sec_GFP_icd; bleach_per_sec_GFP_rpos])
% merge the mcherrys
mean_bleach_per_sec_mCherry_l31=mean([ bleach_per_sec_mCherry_L31_1; bleach_per_sec_mCherry_L31_2])
std_bleach_per_sec_mCherry_l31=std([ bleach_per_sec_mCherry_L31_1; bleach_per_sec_mCherry_L31_2])

%gfp
mean_bleach_per_illum_GFP_rpos=mean(bleach_per_illum_GFP_rpos)
std_bleach_per_illum_GFP_rpos=std(bleach_per_illum_GFP_rpos)
mean_bleach_per_illum_GFP_icd=mean(bleach_per_illum_GFP_icd)
std_bleach_per_illum_GFP_icd=std(bleach_per_illum_GFP_icd)
%merge gfp all together
mean_bleach_per_illum_GFP_merge=mean([ bleach_per_illum_GFP_icd; bleach_per_illum_GFP_rpos])
std_bleach_per_illum_GFP_merge=std([ bleach_per_illum_GFP_icd; bleach_per_illum_GFP_rpos])
% merge the mcherrys
mean_bleach_per_illum_mCherry_l31=mean([ bleach_per_illum_mCherry_L31_1; bleach_per_illum_mCherry_L31_2])
std_bleach_per_illum_mCherry_l31=std([ bleach_per_illum_mCherry_L31_1; bleach_per_illum_mCherry_L31_2])



