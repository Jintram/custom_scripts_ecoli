% *** Extrinsic Noise project ***
% Script to determine whether the measured explained fraction in intrinsic
% noise could also have been observed if the coupling of mu to CFP and YFP
% had been identical 
% YFP and CFP refer to production rate (not concentration)
%
% see also notebook p.43f
% written by NW 2015-01-21
%
% INTIRNSIC NOISE BLUBB USE EXPEIRMENTAL VALUE??? OR DONT DETERMINE FRAC??
%
%
% **********************************
% SIGNIFICANCE TEST
% **********************************
% Null-Hyptothesis: CFP and YFP have same coupling strength to mu
% Alternative Hypothesis: the coupling is different
%
% Experimental observation: The explained fraction of intrinsic noise was
%        XX% (to be specified for each experiment, mostly roughly 3%),
%        how likely is this (or a more extreme) fraction if the 
%        Null-hypothesis (= equal coupling) was true?
%
% ---------------------------
% Practical Implementation:
% ---------------------------
% (1) We don't know how to sample representative distributions and the coupling
% between these distributions: "What is a real distribution supposed to
% look like"
% Make some simplifying assumptions:
% * use the measured corrrelations and extr/intr noise to define a sample dataset
% * Interprete coupling as correlation (not even sure if that's a
%    simplification)
% * Assume linearity (this already is somewhat based on the model. but the
%    real system is close to linear - TODO blubb: acetate is less linear: does
%    it work there also?)
% * Assume gaussian distributions of mu, CFP, YFP
%
% (2) Create a dataset of the same size and correlations and noise strength
%    as the experimental data (details see below)
%
% (3) Calculate explained fraction via KDE (Binning not implemented yet)
%
% (4) Repeat (2)&(3) ca 10.000(?) times to extract the probability to have
%    an expl frac of >XX%


% **************************************************************************
% Summary of experimental data
% **************************************************************************
% 2012-05-08pos2 Maltose full       [time: 2.000sec or 10.000 samples
% numdata=4225;
% noisemu=0.1937;
% noiseyfp=0.4343;
% noisecfp=0.4151; % use min/max/min of these 2 different noises
% noise2_intr=0.0485;
% corryfpmu=0.2817; % use min/max/mean fo the 2 different corrs
% corrcfpmu=0.1715;
% experimental_frac_intr_expl=0.029646;


% 2012-06-15pos4 Acetate full
% numdata=1900;
% noisemu=0.2479;
% noiseyfp=0.5792;
% noisecfp=0.5266; % use min/max/min of these 2 different noises
% noise2_intr=0.0712;
% corryfpmu=0.5295; % use min/max/mean fo the 2 different corrs
% corrcfpmu=0.4966;
% experimental_frac_intr_expl=0.025599;


% 2012-07-26 Maltose basal
%numdata=1947; 
%noisemu=0.1705;  
%noiseyfp=1.2635
%noisecfp=1.1141
%noise2_intr=0.8084; 
%corryfpmu=0.0961;  
%corrcfpmu=0.0631;
%experimental_frac_intr_expl=0.0060775;
% **************************************************************************


% **************************************
% *** ADJUST ***
% **************************************
numrepeats=10000;   % # of times the sampling of vectors & expl frac calc is repeated
SHOWOUTPUT=0;   % display some data about sampling vectors
numdata=1947;   % # of data points the experiment has
noisemu=0.1705;    % noise in growth rate (std/mean)
noiseyfp=mean([1.1141,1.2635]);   % noise in YFP production rate [blubb Note: currently this 
                % refers to total noise - contrary to "all noise which is
                % due to sources other than mu"
                % This could become relevant for different coupling
                % strengths
                % default: choose mean of experim. noise(CFP),noise(YFP)
noisecfp=noiseyfp;
noise2_intr=0.8084;  % squared(!) intrinsic noise. "Taken as a given constant in a certain experiment!!!"

corryfpmu=0.0961;  % correlation btw YFPrate and mu. Choose either min, max or
                % mean of [Corr(YFP,mu),Corr(CFP,mu)]
                % Also test more extreme values
corrcfpmu=0.0631;%corryfpmu;

experimental_frac_intr_expl=0.0060775;   % experimentally measured intr expl fraction (KDE)
% **************************************




% prepare results vectors
noise_intr_expl_abs_vec=nan(numrepeats,1); % expl intr noise (abs! not frac) via KDE.
noise_intr_expl_frac_vec=nan(numrepeats,1); % expl intr frac
% for the fun also the extrinsic ones
noise_extr_expl_abs_vec=nan(numrepeats,1); % expl extr noise (abs! not frac) via KDE.
noise_extr_expl_frac_vec=nan(numrepeats,1); % expl extr frac

% ------------------------------------------------------------------
% loop over all replicates
% ------------------------------------------------------------------
tic
for num=1:numrepeats

% ------------------------------------------------------------------
% create vectors
% ------------------------------------------------------------------
    muvecgauss=randn(numdata,1); % gaussian distribution of mu's
    muvec=1+noisemu*muvecgauss;  % adjust mean (=1) and noise

    yfpvecgauss=corryfpmu*muvecgauss+sqrt(1-corryfpmu^2)*randn(numdata,1); % gaussian for
                                % yfp with defined Corr(YFP,mu)
    yfpvec=1+noiseyfp*yfpvecgauss; % adjust mean and variance

    cfpvecgauss=corrcfpmu*muvecgauss+sqrt(1-corrcfpmu^2)*randn(numdata,1); % gaussian for
                                % cfp with defined Corr(YFP,mu)
    cfpvec=1+noisecfp*cfpvecgauss; % adjust mean and variance


    % some output - deactivate when numrepeats>>1
    if SHOWOUTPUT
        disp(' ')
        disp('-----------------')
        disp(['Constructed vectors:'])
        disp('-----------------')
        disp(['mu: mean=' num2str(mean(muvec)) ',    noise=' num2str(std(muvec)/mean(muvec)) '.    Theo input: noise(mu)=' num2str(noisemu)])
        disp(['YFP: mean=' num2str(mean(yfpvec)) ',    noise=' num2str(std(yfpvec)/mean(yfpvec)) '.    Theo input: noise(yfp)=' num2str(noiseyfp)])
        disp(['CFP: mean=' num2str(mean(cfpvec)) ',    noise=' num2str(std(cfpvec)/mean(cfpvec)) '.    Theo input: noise(cfp)=' num2str(noisecfp)])
        disp('Correlations:')
        disp(['Corr(YFP,mu)=' num2str(corr(yfpvec,muvec)) '.   Theo input: ' num2str(corryfpmu)]);
        disp(['Corr(CFP,mu)=' num2str(corr(cfpvec,muvec)) '.   Theo input: ' num2str(corrcfpmu)]);
        disp(' ')
    end

% ------------------------------------------------------------------
% calculate explained fraction via KDE
% ------------------------------------------------------------------
    % CAREFUL! THIS ARTIFICIAL DATA HAS THE WRONG INTRINSIC AND EXTRINSIC
    % NOISE VALUES! ONLY THE CORRELATIONS (CFP,mu) and (YFP,mu) are
    % CORRECT. THEREFORE ONLY THE noise_XX_expl_abs ARE CORRECT. USE THE
    % EXPERIMENTAL NOISE_EXTR/INTR VALUES
    
    %[noise_extr,noise_extr_expl_abs, noise_extr_expl_frac, noise_intr, noise_intr_expl_abs, ...
    %    noise_intr_expl_frac] = Func_NoiseDecomposition_viaKDE2(cfpvec, yfpvec, muvec);
    [~,noise_extr_expl_abs, ~, ~, noise_intr_expl_abs, ~] =  ...
        Func_NoiseDecomposition_viaKDE2(cfpvec, yfpvec, muvec);
          
    % calculate the explained fractions (assuming cst extr & intr noise)
    noise_intr_expl_frac_calc=noise_intr_expl_abs/noise2_intr;
    noise2_extr_calc=mean([noisecfp^2,noiseyfp^2])-noise2_intr;
    noise_extr_expl_frac_calc=noise_extr_expl_abs/noise2_extr_calc;
    
    % combine data
    noise_intr_expl_abs_vec(num,1)=noise_intr_expl_abs;
    noise_intr_expl_frac_vec(num,1)=noise_intr_expl_frac_calc;
    noise_extr_expl_abs_vec(num,1)=noise_extr_expl_abs;
    noise_extr_expl_frac_vec(num,1)=noise_extr_expl_frac_calc;
    
    
    % some update where we are
    if mod(num,500)==0
        disp(['Finished repeat ' num2str(num)]);
    end
    
end
toc

% ------------------------------------------------------------------
% plots and p value
% ------------------------------------------------------------------

idx=find(isnan(noise_intr_expl_frac_vec));
if ~isempty(idx)
    error('nan values!')
end
num_extreme_events=find(noise_intr_expl_frac_vec>=experimental_frac_intr_expl);
num_extreme_events=length(num_extreme_events);
pvalue=num_extreme_events/numrepeats;

meanintrexplfrac=mean(noise_intr_expl_frac_vec);
meanextrexplfrac=mean(noise_extr_expl_frac_vec);

% intr noise
figure
clf; hold on;
set(gcf,'WindowStyle','docked')
hist(noise_intr_expl_frac_vec,20);
xlabel('Expl frac of \eta_{intr}^2')
ylabel('frequency')
title(['Expl frac of \eta_{intr}^2: mean=' num2str(meanintrexplfrac) '.  p-value for >' ... 
    num2str(experimental_frac_intr_expl) ':  ' num2str(pvalue)])


% extr noise
figure
clf; hold on;
set(gcf,'WindowStyle','docked')
hist(noise_extr_expl_frac_vec,20);
xlabel('Expl frac of \eta_{extr}^2')
ylabel('frequency')
title(['Expl frac of \eta_{extr}^2: mean=' num2str(meanextrexplfrac) ])

