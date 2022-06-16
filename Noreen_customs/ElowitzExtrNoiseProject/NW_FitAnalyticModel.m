% SEVERAL STEPS TO EXTRACT THE FITTING PARAMETERS FOR THE ANALYTICAL MODEL
% (cf notes p. 153ff)

% load all schnitzcells.mat files
%load '\\biofysicasrv\Users2\Walker\ExperimentsCollectedData\schnitzcellsCollectionForNoiseExplByMu20120913.mat';

% create pseudo p-structure (only belangrijk feature: create save
% directory)
p = DJK_initschnitz('plots','2012-xx-yy','e.coli.AMOLF','rootDir','\\biofysicasrv\Users2\Walker\ExtrNoiseProjectAnalyticModel\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','y','fluor2','c','fluor3','none');
% *** TODO *** HOW TO CHOOSE IF SCATTERPLOT AND BRANCHDATA NOISE IS VERY
% DIFFERENT! **** WHEN DOES FITTIME COME IN? -> can destroy the length of a
% lot of branches if set from useForPlot!

% ************************************************
% *****SCATTER CALCULATION ***************************
% schnitzUseName='schnitzcells_malt20120516pos6_basal';  %schnitzcells selection
% %schnitzcell selection might have fight time set -> removes a lot of (too
% %much?) data. maybe load just s_rm
% eval(['schnitzUse=' schnitzUseName ';']);
% fitTime = [0 5000]; fitTime = fitTime + [2 -2];
% **** MAYBE BETTER FOR BRANCH CALCULATIONS ***

schnitzUseName='schnitzcells_malt20120516pos6_basal';
fullName=['rm_' schnitzUseName]; fullTimeName=['fitTime_' schnitzUseName];
eval(['schnitzUse=' fullName ';']);
eval(['fitTime=' fullTimeName ';']);
% *********************************************

% ****
% fields to use for crosscorr resp. crosscov (choose 2 fields, if you want
% autocorr then choose same field twice
field1='dY5_cycCor';%'dY5_sum_dt_s_cycCor'; %'dY5_cycCor'; %'muP15_fitNew_cycCor'
field2='muP15_fitNew_cycCor';
% time field associated to this field (choose rarer time field)
mytimefield='dY5_time';  %'Y_time'; %'dY5_time'; %'time_atdY';
% ****

% *** decide how to make average!!!  *** all data (s_rm) or only within
% fitTime (more work)
% normalize used fields by their average (otherwise, crosscovariance is
% weird. it also fits better to the model).
allfield1=[]; allfield2=[];
%mytime=[]; %debug
for i=1:length(schnitzUse)
    if schnitzUse(i).useForPlot==1
        % ** only take data into account that is within fitTime **
        idx1=min(find(schnitzUse(i).(mytimefield)>=fitTime(1))); %fT
        idx2=max(find(schnitzUse(i).(mytimefield)<=fitTime(2))); %fT
        if ~isempty(idx1) & ~isempty(idx2) & (idx1<=idx2) %fT
            allfield1=[allfield1 schnitzUse(i).(field1)(idx1:idx2)];
            allfield2=[allfield2 schnitzUse(i).(field2)(idx1:idx2)];
            % mytime=[mytime schnitzUse(i).(mytimefield)(idx1:idx2)]; %debug
        end %fT
    end
end            
%figure; hist(mytime,30) %debug
%allfield1=[schnitzUse.(field1)]; % does not check useForPlot
notnanidx=find(~isnan(allfield1));
allfield1=allfield1(notnanidx);
meanfield1=mean(allfield1);
%allfield2=[schnitzUse.(field2)]; % does not check useForPlot
notnanidx=find(~isnan(allfield2));
allfield2=allfield2(notnanidx);
meanfield2=mean(allfield2);

schnitzUseNorm=schnitzUse; % ..Norm: averages set to 1. If not wanted, continue with SchnitzUse
for i=1:length(schnitzUseNorm)
    % on right hand side of equaiton: use SchnitzUse (not ..Norm!!).
    % Reason: If autocorrelation is calculated, field1=field2 is otherwise
    % normalized twice! (alternative: insert check if field1==field2)
    schnitzUseNorm(i).(field1)=schnitzUse(i).(field1)/meanfield1;
    schnitzUseNorm(i).(field2)=schnitzUse(i).(field2)/meanfield2;
   
end 

% just a test *** start ***
allfield1norm=[]; allfield2norm=[];
for i=1:length(schnitzUseNorm) % ignores fitTime!!!
    if schnitzUseNorm(i).useForPlot==1
        allfield1norm=[allfield1norm schnitzUseNorm(i).(field1)];
        allfield2norm=[allfield2norm schnitzUseNorm(i).(field2)];
    end
end    
% allfield1=[schnitzUseNorm.(field1)];  % does not check useForPlot
notnanidx=find(~isnan(allfield1norm));
allfield1norm=allfield1norm(notnanidx);
meanfield1Norm=mean(allfield1norm);
% allfield2=[schnitzUseNorm.(field2)];  % does not check useForPlot
notnanidx=find(~isnan(allfield2norm));
allfield2norm=allfield2norm(notnanidx);
meanfield2Norm=mean(allfield2norm);
% *** end ***

if strcmp(mytimefield,'Y_time')==0 % do not add branchdata of redundant fields ->errormessage
    if strcmp(field1,field2)==0
        branchData = DJK_getBranches(p,schnitzUseNorm,'dataFields',{ mytimefield 'Y_time'  field1 field2 }, 'fitTime', fitTime); 
        branches = DJK_addToBranches_noise(p, branchData,'dataFields',{ mytimefield 'Y_time' field1 field2 });
    else
        branchData = DJK_getBranches(p,schnitzUseNorm,'dataFields',{ mytimefield 'Y_time'  field1 }, 'fitTime', fitTime); 
        branches = DJK_addToBranches_noise(p, branchData,'dataFields',{ mytimefield 'Y_time' field1  });
    end
else
    if strcmp(field1,field2)==0
        branchData = DJK_getBranches(p,schnitzUseNorm,'dataFields',{ mytimefield   field1 field2 }, 'fitTime', fitTime); 
        branches = DJK_addToBranches_noise(p, branchData,'dataFields',{ mytimefield  field1 field2 });
    else
        branchData = DJK_getBranches(p,schnitzUseNorm,'dataFields',{ mytimefield   field1 }, 'fitTime', fitTime); 
        branches = DJK_addToBranches_noise(p, branchData,'dataFields',{ mytimefield  field1  });
    end
end        
trimmed_branches = DJK_trim_branch_data(branches);
branch_groups = DJK_divide_branch_data(trimmed_branches);
myname = schnitzUseName;
ccData=NW_extractCrossCorrData(p, branch_groups, ['noise_' field1], ['noise_' field2],'selectionName',myname, 'showPlots',0);
% for better postprocessing, change names:
T=ccData(:,1);
cc=ccData(:,2);
ccN=ccData(:,3);
ccW=1./ccN;  % weights are inverted to noise
cv=ccData(:,4);
cvN=ccData(:,5);
cvW=1./cvN;

idx=find(T==0);
disp('  ')
disp(['tau=0:  crosscov=' num2str(cv(idx))  '    sqrt(crosscov)=noise: ' num2str(sqrt(cv(idx)))]);
disp(['tau=0:  crosscorr=' num2str(cc(idx))]);

if strcmp(field1,field2)==1
    % compare with result of direct covariance without weighing+branchstructure
    % **** ignores fitTime ***
    scattercov=cov(allfield1norm, allfield2norm); %(1,2) resp. (2,1) is covariance.  (1,1) is variance of field1
    disp(' ')
    disp(['covariance from scatter plot (no weighing)=' num2str(scattercov(1,2))]);
    disp(['noise=sqrt(cov) (normed!) from scatter plot=' num2str(sqrt(scattercov(1,2)))]);
    disp(['If contradictory to histograms, see comments in code.']);
    % ************* NOTE!! ************
    % if the data from scatter plots gives totally different data compared to
    % the histograms (/schnitzcells/) be aware that:
    % - for histograms, only data points within the fitTime are used
    % - here, all datapoints (s_rm) or all cells which are BORN within fitTime
    % (s_rm_fitTime) are used -> more data, especially noisy late data.
    % if the restriction of fitTime is lifted for the histograms, the result is
    % identical to s_rm resp. s_rmfitTime (NW 2012-10)
    % *********************************
end

%% **************************************************************************
% SET SEVERAL OF ANALYTIC MODEL PARAMETERS
% CALCULATE OTHER VARIABLES BY SOLVING USING COVARIANCE AT TAU=0
% (Nomenclature: p.157)
%
% Variables to calculate:
% T_cG, T_yG, Theta_G, Theta_E, Theta_mu, Theta_cy (Theta_c and Theta_y are set equal. m stands for mu)
clear T_cG T_yG Theta_G Theta_E Theta_mu Theta_cy
%
% **** Values extracted from fits ****
beta_G=1.14;  % decay time global noise
beta_mu=1.14; % decay mu specific noise. maybe use decay time according to frameFitRange! (15frames)
R_mc=0.073;   %lin fit values (dY_cycCor etc)
R_my=0.078;
R_cy=0.57;
R_mm=0.10;
R_cc=1.38;
R_yy=1.86;

% maltose20120516pos6 basal
%beta_G=1.14;  % decay time global noise
%beta_mu=1.14; % decay mu specific noise. maybe use decay time according to frameFitRange! (15frames)
%R_mc=0.073;   %lin fit values (dY_cycCor etc)
%R_my=0.078;
%R_cy=0.57;
%R_mm=0.10;
%R_cc=1.38;
%R_yy=1.86;

% maltose20120726pos4 basal
%beta_G=0.56;  % decay time global noise
%beta_mu=0.56; % decay mu specific noise. maybe use decay time according to frameFitRange! (15frames)
%R_mc=0.016;   %lin fit values (dY_cycCor etc)
%R_my=0.027;
%R_cy=0.74;
%R_mm=0.053;
%R_cc=1.27;
%R_yy=1.71;

% acetate20120620pos2
%beta_G=0.24;  % decay time global noise
%beta_mu=0.24; % decay mu specific noise
%R_mc=0.047;   %lin fit values (dY_cycCor etc)
%R_my=0.057;
%R_cy=0.18;
%R_mm=0.054;
%R_cc=0.19;
%R_yy=0.25;

% **** Values trivially set to certain value ****
T_mG=1;      % Transmission global->mu
T_cE=1;      % Transmission protein->c

% **** Further assumptions ****
beta_E=beta_G; %or 60/9 -> explanation? measurement errors? % decay time protein noise
T_yE=T_cE;      % transmission protein->y (strong assumption!)
beta_c=60/9;    % 1/decay time intrinsic noise cfp (literature. check value. rosenfeld2005?)
beta_y=beta_c;  % 1/decay time intrinsic noise yfp (literature)
% Theta_c=Theta_y -> "Theta_cy"  same intrinsic noise strength for yfp and cfp (strong assumption)

% **** simplify parameters by merging the (cf p.157) ***
a=T_mG/(2*beta_G);
bG=1/(2*beta_G);
c=T_cE*T_yE/(2*beta_E); % assumes T_cE=T_yE!
d=T_mG*T_mG/(2*beta_G);
bMu=1/(2*beta_mu);
bCY=1/(2*beta_c); % assumes beta_c=beta_y!
% **** equations ****
% (1) R_mc=a*T_cG*(Theta_G)²
% (2) R_my=a*T_cG*(Theta_G)²
% (3) R_cy=T_cG*T_T_yG*(Theta_G)²*bG
% (4) R_mm=d*(Theta_G)²+(Theta_mu)²*bMu
% (5) R_cc=(T_cG)²*(Theta_G)²*bG + c*(Theta_E)² + bCY*(Theta_cy)²
% (6) R_yy=(T_yG)²*(Theta_G)²*bG + c*(Theta_E)² + bCY*(Theta_cy)²

% **** Calculate ****
klammer=(R_mc/R_my)^2-1;
T_yG=(R_cc-R_yy)*a/(klammer*R_my*bG);
T_cG=R_mc/R_my*T_yG;
Theta_G=sqrt(R_mc/(a*T_cG));
differenz=R_mm-d*(Theta_G^2);
Theta_mu=sqrt(differenz/bMu); clear differenz;
differenz=R_cy-T_cG*T_yG*(Theta_G^2)*bG;
Theta_E=sqrt(differenz/c); clear differenz;
differenz=R_cc-(T_cG^2)*(Theta_G^2)*bG-c*(Theta_E^2);
Theta_cy=sqrt(differenz/bCY); clear differenz;

% *** Convert (new convention) transmissions and noise from global source
% (G) such that T_cG=1 (better comparability to 'E')
newTheta_G=Theta_G*T_cG;
T_cG=1;
T_yG=T_yG*Theta_G/newTheta_G;
T_mG=T_mG*Theta_G/newTheta_G;
Theta_G=newTheta_G;
% redo for plotting abbreviations
a=T_mG/(2*beta_G);
bG=1/(2*beta_G);
c=T_cE*T_yE/(2*beta_E); % assumes T_cE=T_yE!
d=T_mG*T_mG/(2*beta_G);
bMu=1/(2*beta_mu);
bCY=1/(2*beta_c); % assumes beta_c=beta_y!
% ***

% ***** some output ****
disp(' ')
disp('Parameter Values')
disp('-----------------------')
disp(' *** Trivial: *** ')
disp(['T_mG=' num2str(T_mG)])
disp(['T_cE=' num2str(T_cE)])

disp(' *** Assumed *** ')
disp('beta_E=beta_G')
disp('T_yE=T_cE')
disp('beta_c=9/60')
disp('beta_y=beta_c')

disp(' *** Extracted From Plots: *** ')
disp(['beta_G=' num2str(beta_G)])
disp(['beta_mu=' num2str(beta_mu)])
disp(['R_mc=' num2str(R_mc)])
disp(['R_my=' num2str(R_my)])
disp(['R_cy=' num2str(R_cy)])
disp(['R_mm=' num2str(R_mm)])
disp(['R_cc=' num2str(R_cc)])
disp(['R_yy=' num2str(R_yy)])

disp(' *** Calculated: *** ')
disp(['T_yG=' num2str(T_yG)])
disp(['T_cG=' num2str(T_cG)])
disp(['Theta_G=' num2str(Theta_G)])
disp(['Theta_mu=' num2str(Theta_mu)])
disp(['Theta_E=' num2str(Theta_E)])
disp(['Theta_cy=' num2str(Theta_cy)])




%plot crosscorrelations
% ****** if not identical to CC aus branchplot (DJK_standard_eror_store...)
% -> there, for each(!) branch, normalization is performed with
% sqrt(Var(x)*Var(y)) -> can be different than averaging Var first and then
% dividing the crosscov
% ************************
T=-10:0.1:10; T0=0:0.1:10;
R_mcTau=a*T_cG*(Theta_G^2)*exp(-beta_G*abs(T))/sqrt(R_mm*R_cc);  %including normalization with sigma at tau=0
R_myTau=a*T_yG*(Theta_G^2)*exp(-beta_G*abs(T))/sqrt(R_mm*R_yy);  %including normalization with sigma at tau=0
R_cyTau=(T_cG*T_yG*(Theta_G^2)*bG*exp(-beta_G*abs(T)) + c*(Theta_E^2)*exp(-beta_E*abs(T))) / sqrt(R_cc*R_yy);
R_mmTau=(d*(Theta_G^2)*exp(-beta_G*T0)+(Theta_mu^2)*bMu*exp(-beta_mu*T0) ) / (R_mm);
R_ccTau=(T_cG^2*Theta_G^2*bG*exp(-beta_G*T0)+c*Theta_E^2*exp(-beta_E*T0)+bCY*Theta_cy^2*exp(-beta_c*T0) ) / R_cc;
R_yyTau=(T_yG^2*Theta_G^2*bG*exp(-beta_G*T0)+c*Theta_E^2*exp(-beta_E*T0)+bCY*Theta_cy^2*exp(-beta_y*T0) ) / R_yy;


h=figure(6);
clf
grid on
hold on
plot(T,R_mcTau,'b')
plot(T,R_myTau,'r')
plot(T,R_cyTau,'Color',[0.5 0.5 0.5]) %gray
plot(T0,R_mmTau,'Color',[0 0.8 0.8])  %dark cyan
plot(T0,R_ccTau,'Color',[0.9 0 0.8])  %dark magenta
plot(T0,R_yyTau,'Color',[1 0.8 0])    %orange

legend('R\_mc','R\_my','R\_cy','R\_mm','R\_cc','R\_yy')

