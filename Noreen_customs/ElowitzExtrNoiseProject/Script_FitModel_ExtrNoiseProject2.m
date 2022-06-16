% Extrinsic noise project: Perform the model fitting
%
% Version: 2015-02-26
% For detail on the model see handwritten notebook p. 59ff
%
% FITTED PARAMETERS: a_YG,a_CG,a_muG: coupling constants to global noise
%                                   we set a_muG=1 (arbitrary
%                                   normalization)
%                   a_YP, a_CP: coupling constants to protein noise
%                   a_XX (=a_YY=a_CC): coupling constant to intrinsic noise
%                   var_S:  cst mu-specific noise (meaning: variance)
%                                   "offset" (var_S=\eta_S^2)
%                   var_P0: cst part protein specific noise
%                   var_C0, var_Y0 : cst part YFP & CFP specific intr noise
%
%
% DATA:
% Run the first cell in
% "\\biofysicasrv\Users2\Walker\ExtrinsicNoise\Analysis\Plot_ExplFrac_Covs_vs_NoiseGrowth\
%       script_plot_explfractions_etc_vs_noisegrowth.m "
% Creates all vectors: noise, mu, description, etc
% Be careful with changing any entries! (They are fixed assigned values,
% not loaded automatically!)
%
% For fitting use the subset of datapoints: idxfullcst_unique
%       = all fully induced + constitutive promoter datasets, which are
%       unique (=excludes the 2nd noisier point of 2012-05-24 with bimodal growth distr)

% **********************************************
% For formulas and enumeration check handwritten notes p. 59ff
%   a_muG=1. Var_G is the varying noise source (and Var_P etc scale with
%   Var_P ..., note that Var_P0 etc doesn't scale)
%   Var_mu = var_G + var_S
%
%   (I) Var_C = (a_CG^2+a_CP^2+a_XX^2)*Var_G + Var_P0 + Var_C0
%   (II) Cov(C,mu) = ....
%   ....
% **********************************************

clear a_YG a_CG a_YP a_CP a_XX var_S var_P0 var_Y0 var_C0

% get data together in easily handable way
myidx=idxfullconst_unique;

varmuforFit=noisemu(myidx).^2;
varCforFit=noisedC(myidx).^2;
varYforFit=noisedY(myidx).^2;
cov_CYforFit=cov_dCdY(myidx);
cov_CmuforFit=cov_dCmu(myidx);
cov_YmuforFit=cov_dYmu(myidx);


% ------------------------------------
% (1) & (2) preliminary fits of Cov(C,mu), Cov(Y,mu)
% ------------------------------------
% first Cov(dY,mu) fit
yfitP=polyfit(varmuforFit,cov_YmuforFit,1);  % P=prelim
yslopeP=yfitP(1);   %1.2722
yoffsetP=yfitP(2);
var_S_viaYfit=-yoffsetP/yslopeP;

% first Cov(dC,mu) fit
cfitP=polyfit(varmuforFit,cov_CmuforFit,1);  % P=prelim
cslopeP=cfitP(1);   %0.9327
coffsetP=cfitP(2);
var_S_viaCfit=-coffsetP/cslopeP;

% ------------------------------------
% (3) & (4) average the Var_S fits and refit the Cov(dY,mu), Cov(dC,mu) slopes
% VARIABLES: a_CG, a_YG, var_S
% ------------------------------------
% average the var_S-fits
var_S_fit=mean([var_S_viaCfit,var_S_viaYfit]);

% refit with only slope as free parameter.
% largely copied form curve fitting toolbox
ft_ = fittype({['(x-' num2str(var_S_fit) ')']}, 'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a'});
% Fit this model using new data [skipped the inf/nan test]
cf_y = fit(varmuforFit(),cov_YmuforFit(),ft_);
yslope=cf_y.a;      %1.2678 (very similar to above! :-) )
cf_c = fit(varmuforFit(),cov_CmuforFit(),ft_);
cslope=cf_c.a;      %0.9359 (very similar to above! :-) )


% Rename
a_YG=yslope;
a_CG=cslope;
var_S=var_S_fit;


% ------------------------------------
% (5) fit all slopes.
% name all the slopes and offsets in a simple way.
% use Var(C) and Var(Y) and Cov(Y,C)
% VARIABLES: a_CP
% ------------------------------------
% (I) Var(C)
dummy=polyfit(varmuforFit,varCforFit,1);  
D=dummy(1); D0=dummy(2); clear dummy
% (II) Cov(C,mu)
E=a_CG; % already done (could not be done at beginning because joint var_S needed)
% (III) Var(Y)
dummy=polyfit(varmuforFit,varYforFit,1);  
F=dummy(1); F0=dummy(2); clear dummy
% (IV) Cov(Y,mu)
G=a_YG; % already done
% (V) Cov(Y,C)
dummy=polyfit(varmuforFit,cov_CYforFit,1);  
H=dummy(1); H0=dummy(2); clear dummy

% calculate a_CP
firstterm=0.5*(-(E*E-G*G+F-D));
secondterm=0.5*sqrt((E*E-G*G+F-D)^2+4*(H-G*E)^2);
% only one root of a_CP=firstterm+-secondterm should work out
if secondterm<firstterm
    error('Both a_CP roots are positive!!!')
else
    a_CP=sqrt(firstterm+secondterm);
end


% ------------------------------------
% (6) use Var(C) and Var(Y)
% VARIABLES: a_XX, a_YP
% ------------------------------------
a_XX=sqrt(D-a_CG^2-a_CP^2);
a_YP=sqrt(F-a_XX^2-a_YG^2);


% ------------------------------------
% (7) name some coefficents combinations in a simpler way.
% now calculate offsets.
% uses Cov(Y,C)
% VARIABLES: var_P0
% ------------------------------------
gamma=a_YG*a_CG+a_YP*a_CP;
delta=a_CG^2+a_CP^2+a_XX^2;
phi=a_YG^2+a_YP^2+a_XX^2;

var_P0=H0+gamma*var_S;

% ------------------------------------
% (8) uses Var(C)
% VARIABLES: var_C0
% ------------------------------------
var_C0=D0+delta*var_S-var_P0;

% ------------------------------------
% (9) uses Var(Y)
% VARIABLES: var_Y0
% ------------------------------------
var_Y0=F0+phi*var_S-var_P0;

% ----------------------------------------------
% FINISHED FITTING
% ----------------------------------------------

% Prepare some lines for the plotting

%range of var_mu plotted
mumin=max(min(varmuforFit)-0.01,var_S);
Var_mu_forPlot=mumin:0.001:max(varmuforFit)+0.01;

% Var & Cov -> Linear linees used for fitting parameters
Var_C_Fitted=delta*(Var_mu_forPlot-var_S)+var_P0+var_C0;
Cov_Cmu_Fitted=a_CG*(Var_mu_forPlot-var_S);
Var_Y_Fitted=phi*(Var_mu_forPlot-var_S)+var_P0+var_Y0;
Cov_Ymu_Fitted=a_YG*(Var_mu_forPlot-var_S);
Cov_YC_Fitted=gamma*(Var_mu_forPlot-var_S)+var_P0;

% Corrs
Corr_Cmu_Fitted=Cov_Cmu_Fitted./(sqrt(Var_mu_forPlot).*sqrt(Var_C_Fitted));
Corr_Ymu_Fitted=Cov_Ymu_Fitted./(sqrt(Var_mu_forPlot).*sqrt(Var_Y_Fitted));
Corr_YC_Fitted=Cov_YC_Fitted./(sqrt(Var_Y_Fitted).*sqrt(Var_C_Fitted));

% Extr & Intr squared noise
Extrnoise2_Fitted=Cov_YC_Fitted;
Intrnoise2_Fitted=0.5*((a_YP-a_CP)^2+(a_YG-a_CG)^2+2*a_XX^2)*(Var_mu_forPlot-var_S) ...
    +0.5*var_C0+0.5*var_Y0;

% Explained noise & Explained fraction
Extrnoise2Explained_Fitted=a_CG*a_YG*(Var_mu_forPlot-var_S).^2./Var_mu_forPlot;
Extrnoise2FracExpl_Fitted=Extrnoise2Explained_Fitted./Extrnoise2_Fitted;
Intrnoise2Explained_Fitted=0.5*(a_CG-a_YG)^2*(Var_mu_forPlot-var_S).^2./Var_mu_forPlot;
Intrnoise2FracExpl_Fitted=Intrnoise2Explained_Fitted./Intrnoise2_Fitted;

% -----------------------------
% FOR PLOTTING: PLOT THESE LINES INTO Script_ModelFit....
% -----------------------------

%%
% ---------------------------------------
% fit the basal covariance cov(Y,C)
% ---------------------------------------

varmuforFitBasal=noisemu(idxbasal_unique).^2;
cov_CYforFitBasal=cov_dCdY(idxbasal_unique);

fitBasal=polyfit(varmuforFitBasal,cov_CYforFitBasal,1);
slopebasal=fitBasal(1);  
offsetbasal=fitBasal(2);

Cov_YC_FittedBasal=slopebasal*Var_mu_forPlot+offsetbasal;
%the offset is not converted into the noise parameters var_P0 etc


%%

% some plotting

%figure
clf 
hold on
plot(varmuforFit,cov_CdYforFit,'.','MarkerSize',20)
plot(xx,aY*aC*(xx-nS)+nYC,'k')

