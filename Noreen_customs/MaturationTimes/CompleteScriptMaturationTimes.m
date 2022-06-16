

% Run this script to perform all steps for maturation time calculation
% along branches

% schnitzcell software must have been run before. 
% necessary data in workspace: p, s_rm (schnitzcells with removed sick schnitzes)



%% *****************************************************
% CREATE NEW TIMEFIELD - ONLY RUN ONCE !!!

% **** ADJUST ***************************
% %timeshift wrt 1st phase image of pos7
% %for pos7 : =0    % 2013-06-25
% %for pos3 : =+15

% timeshiftinittime=138/60;

% *************************************

%for i=1:length(schnitzcells)
%    schnitzcells(i).G_time_wrtPos2phase=schnitzcells(i).G_time+timeshiftinittime;
%end
%NW_saveSchnitzcells(p,schnitzcells);


% % RUN SCHNITZCELLS RESTRICTIONS -> s_rm




%% *****************************************************
% GET LINEAGE DATA
% *** ADJUST **********
timefield='G_time_wrtPos7phase';
%yfield='G5_sum';
yfield='R5_sum';
nrbranches=4;
myp=pos7p;
%myp=pos4p;
myschnitzcells=pos7s_rm;
%myschnitzcells=pos4schnitzcells;
% *****************

[branches, timedata, ydataall] = NW_plotLineagesTimeTrace (myp, myschnitzcells,timefield, yfield, ...
    'extensive',0,nrbranches);

% % remove a bad branch
%badbranch=2;
%branches=branches([1:badbranch-1, badbranch+1:end]);
%ydataall=ydataall([1:badbranch-1, badbranch+1:end],:);


% average fluo curve over time
meanydata=mean(ydataall);

%% *****************************************************
% CALCULATE BLEACHING. CORRECT FOR BLEACHING (ydataallcorr) AND SMOOTH DATA
% (ydataallcorrsmooth)
%
% 
% *** ADJUST ***************************************
% use late time points when bleaching is the only effect (no maturation, no
%                                                    production)
%fittimebleaching=[155 200];
fittimebleaching=[400 450];
% *****************************************************

idxtime=find(timedata>fittimebleaching(1) & timedata<fittimebleaching(2));
timedatasub=timedata(idxtime);
meanydatasub=meanydata(idxtime);

% fit an exponential bleaching function
% (basically copied form "GenerateCode" of Fitting Toolbox)
% check for Nan and Inf values
ok2 = isfinite(timedatasub) & isfinite(meanydatasub);
if ~all( ok2 )
    warning( 'CompleteScriptMaturationTimes: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
end
% suitable start paramters for [a b]
st2= [ 400000 -0.01];
ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
% Fit this model using new data
cfbleach = fit(timedatasub(ok2)',meanydatasub(ok2)',ft2,'Startpoint',st2)

% get fitted line
a=cfbleach.a; bleachrate=-cfbleach.b;
timeaxisbleach=fittimebleaching(1):1:fittimebleaching(end);
yaxisfitbleach=a*exp(-bleachrate*timeaxisbleach);
figure(3);
clf
hold on
set(gcf,'WindowStyle','docked')
plot(timedata,meanydata,'.','MarkerSize',15)
plot(timeaxisbleach,yaxisfitbleach,'r','LineWidth',2)
xlabel('time')
ylabel('total fluo')
title('bleach fit of average fluo data')

% convert bleachrate into bleach factor per illumination event
% !!! ASSUMPTION !!! Time points (only) in fit range are equally spaced!
% (NW 2013-07-10: I don't think that equally spaced actually plays a role)
% Each moment an image is acquired, the cells are bleached by the same
% factor (!!)
% y2 = y1 * exp(-bleachrate * timestep)
%    = y1 * bleachfactor
meantimediff=mean(diff(timedatasub));
bleachfactor=exp(-bleachrate*meantimediff);
disp(' ')
disp(['bleaching rate: '  num2str(bleachrate) '.  bleachfactor per image: ' num2str(bleachfactor)])


% array with all bleaching correction factors (i'th image: 1/bleachfactor^2)
bleachcorrfactors=zeros(size(timedata));
for i=1:length(bleachcorrfactors)
    bleachcorrfactors(i)=1/(bleachfactor^i);
end

% correct for bleaching and smooth data
ydataallcorr=zeros(size(ydataall));
ydataallcorrsmooth=zeros(size(ydataall));
for i=1:size(ydataall,1)
    ydataallcorr(i,:)=ydataall(i,:).*bleachcorrfactors;
    ydataallcorrsmooth(i,:)=smooth(timedata,ydataallcorr(i,:),'lowess');
end
figure(4);
clf
hold on
set(gcf,'WindowStyle','docked')
plot(timedata,ydataall,'.')
title('total fluo. each lineage')
xlabel('time [min]')
ylabel('total fluo')
figure(5);
clf
hold on
set(gcf,'WindowStyle','docked')
plot(timedata,ydataallcorr,'.')
title('total fluo. bleaching corr. each lineage')
xlabel('time [min]')
ylabel('total fluo')
figure(6);
clf
hold on
set(gcf,'WindowStyle','docked')
plot(timedata,ydataallcorrsmooth,'.')
title('total fluo. bleaching corr. smoothed. each lineage')
xlabel('time [min]')
ylabel('total fluo')

%% **********************************************************
% BLEACHING CORRECTION PER LINEAGE

% first smoothed
ydataallsmooth=zeros(size(ydataall));
for i=1:size(ydataall,1)
    ydataallsmooth(i,:)=smooth(timedata,ydataall(i,:),'lowess');
end
ydataallsmoothindivcorr=zeros(size(ydataall));

% *** ADJUST ***************************************
% use late time points when bleaching is the only effect (no maturation, no
%                                                    production)
%fittimebleaching=[160 200];
fittimebleaching=[400 450];
% *****************************************************

idxtime=find(timedata>fittimebleaching(1) & timedata<fittimebleaching(2));
timedatasub=timedata(idxtime);

% loop over all lineages
clear lin
for lin=1:size(ydataallsmooth,1)
    linydatasub=ydataallsmooth(lin,idxtime);

    % fit an exponential bleaching function
    % (basically copied form "GenerateCode" of Fitting Toolbox)
    % check for Nan and Inf values
    ok2 = isfinite(timedatasub) & isfinite(linydatasub);
    if ~all( ok2 )
        warning( 'CompleteScriptMaturationTimes: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
    end
    % suitable start paramters for [a b]
    st2= [ 400000 -0.01];
    ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
    % Fit this model using new data
    cfbleach = fit(timedatasub(ok2)',linydatasub(ok2)',ft2,'Startpoint',st2)

    % get fitted line
    a=cfbleach.a; bleachrate=-cfbleach.b;
    
    % convert bleachrate into bleach factor per illumination event
    % !!! ASSUMPTION !!! Time points (only) in fit range are equally spaced!
    % Each moment an image is acquired, the cells are bleached by the same
    % factor (!!)
    % y2 = y1 * exp(-bleachrate * timestep)
    %    = y1 * bleachfactor
    meantimediff=mean(diff(timedatasub));
    bleachfactor=exp(-bleachrate*meantimediff);
    disp(' ')
    disp(['lineage '  num2str(lin) '. bleaching rate: '  num2str(bleachrate) '.  bleachfactor per image: ' num2str(bleachfactor)])


    % array with all bleaching correction factors (i'th image: 1/bleachfactor^2)
    bleachcorrfactors=zeros(size(timedata));
    for i=1:length(bleachcorrfactors)
        bleachcorrfactors(i)=1/(bleachfactor^i);
    end

    % correct for bleaching and smooth data
    ydataallsmoothindivcorr(lin,:)=ydataallsmooth(lin,:).*bleachcorrfactors;

end
figure(7);
clf
hold on
set(gcf,'WindowStyle','docked')
plot(timedata,ydataallsmooth,'.')
title('total fluo. smooth. each lineage')
xlabel('time [min]')
ylabel('total fluo')
figure(8);
clf
hold on
set(gcf,'WindowStyle','docked')
plot(timedata,ydataallsmoothindivcorr,'.')
title('total fluo. smooth bleaching corr per lineage.')
xlabel('time [min]')
ylabel('total fluo')


%% ********************************************************
% FIT MATURATION TIME KINETICS TO EACH LINEAGE


% *** ADJUST *********** (check FitMaturationTime.m for more info)
% *********
%mytimeOnsetAB=140;
%mytimeMaxFit=160;
%mytimerangeAllFluoMature=[153 154]; %1datapoint vs 2 datapoints [151 154];
        % typically highest point before bleaching dip is chosen
showplot=1;  % 0 or 1    % if restrictToBranches is not empty, this variable will be overwritten
% **************************************************************************
 mytimeOnsetAB=300;   % nice fit pos7 . 25.6.2013 mcherry
 mytimeMaxFit=400;
 mytimerangeAllFluoMature=[380 430];

% mattime is very sensitive to ytimeOnsetAB
% mattime is not sensitive to mytimeMaxFit
% mattime is very sensitive to mytimerangeAllFluoMature(1)  (min. time for av)
%                   large time -> larger mattime
% mattime is a bit dependent on mytimerange...(2)   (max. time for av)
%       tendency: time increase -> mattimeincrease. but not super strict

% ydatause=ydataallsmoothindivcorr;  % indiv. bleaching correction for each branch
 ydatause=ydataallcorrsmooth;  % one common bleaching correction for each

% branch
% correction of individual lineages for bleaching did not really help

% if not empty: only these branches are plotted. all branches are
% calculated
%restrictToBranches=[];%  9.7.pos4: [1 2 3 4 5  8 10 11 14];
restrictToBranches=[1 2 3 5 6]
% **************************************

lineageindices=1:size(ydatause,1); % nr of each different lineage
mattimevec=[];   % maturation time for every lineage
matratevec=[];   % maturation rate for every lineage
matrate_confint_vec=zeros(0,2);  % confidence interval of maturation rates
Fluoimmatvec=[];  % immaturate Fluorophores at start point of fitting for each lineage
Fluoimmat_confint_vec=zeros(0,2);  % confidence interval
Fluototalvec=[];  % total fluorescence (at late-time-plateau) for each lineage)

for lin=1:size(ydatause,1)
    myxdata=timedata;
    myydata=ydatause(lin,:);
    
    if ~isempty(restrictToBranches)
        if ismember(lin,restrictToBranches)
            showplot=1;
        else showplot=0;
        end
    end
  
    %fit
    [matrate,Fimmat,matrate_confint,Fimmat_confint,mattime,Ftotal]=FitMaturationTimes( ...
        lin, myxdata, myydata, mytimeOnsetAB, mytimeMaxFit, mytimerangeAllFluoMature,showplot);
    
    % get data together
    mattimevec=[mattimevec; mattime];
    matratevec=[matratevec; matrate];
    matrate_confint_vec=[matrate_confint_vec; matrate_confint'];
    Fluoimmatvec=[Fluoimmatvec; Fimmat];
    Fluoimmat_confint_vec=[Fluoimmat_confint_vec; Fimmat_confint'];
    Fluototalvec=[Fluototalvec;Ftotal];
end
    
% **** ADJUST *****
% which lineages look nice (good fit, leveling off)
%useLineages=[1 2 3 4 5 8 10 11 14]; %[1 2 3 4];
useLineages=[1 2 3 5 6];
% *************************

% mattimevec;
mattimevec(useLineages)

disp(['mean(mattime)=' num2str(mean(mattimevec(useLineages))) ]);
disp(['stddev(mattime)=' num2str(std(mattimevec(useLineages))) ]);
figure(9)
clf
hold on
hist(mattimevec(useLineages));
set(gcf,'WindowStyle','docked')

