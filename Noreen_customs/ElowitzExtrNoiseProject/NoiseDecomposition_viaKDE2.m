% This script alculates to which part mu-related fluctuations are responsible for
% extrinsic and intrinsic noise
% Total extrinsic noise resp. intrinsic noise is split up according to law
% of total covariance resp. variance
% extr noise: Cov(C,Y|mu)= Cov(<C|mu>,<Y|mu>) + <Cov(C,Y|mu>) = explained + unexplained
% intr noise: Var(C-Y)   = Var(<C-Y|mu>)      + <Var(C-Y|mu)> = explained + unexplained
%
% It's based on KERNEL DENSITY ESTIMATE (core function:
% NW_Probabilities_Grid_viaKDE2.m) get probability distributions.
%
% An alternative method for noise decomposition is binning (use as
% validation of this KDE method to show "robustness of explained fraction
% to chosen method"). script: ExtrIntrNoise_ExplainedByMu_VectorInput.m
% ((de)activate TEXTOUTPUT directly in that function)

% Input vectors are normalized to mean=1 -> if original mean=0 this can
% lead to artefacts!

% ******************
% NO HISTORY DEPENDENCE IMPLEMENTED. FOR ONE TIME STEP IN HISTORY, 
% use ExtrIntrNoise_ExplainedByMu_1TimeptHist.m
% only implemented for Binning method
% ******************

% This script is rather simple and only allows for one parameter set per
% evaluation - to check parameter robustness run this script several times

% ***********************************************************************
% INTIRNSIC NOISE IS NOT CORRECTLY NORMALZIED YET BY THE FACTOR 2 (!!!!)
% ***********************************************************************


%% 0) FACULTATIVE: LOAD A .MAT FILE CONTAINING SCHNITZCELLS AND EXTRACTED VECTORS
% But schnitzname has to be set!
% ************************************************************************
% 
clear myschnitzname
 myschnitzname='R';  % Adjust
% myschnitzname='Acetate full 2012-06-15 pos4';  % Adjust
 load(['\\biofysicasrv\Users2\Walker\ExtrinsicNoise\Data_Collection\Data\' myschnitzname '.mat']);
%

% ************************************************************************
%% 1) SET PLOT/BINNING OPTIONS AND GET DATA TOGETHER
% ************************************************************************

% ----------------------------------------------------------
% ADJUST
% ----------------------------------------------------------
% Plotting options Adjust to KDE output
PLOTMUDISTRIBUTION=0;    % mu (field3) probability distribution
PLOTCONDVALUESEXTR=1;   % extr noise: <C|mu> vs mu plots (same for Y and Cov)
PLOTCONDVALUESINTR=1;   % intr noise <C-Y|mu>
PLOTPROBHEATMAPS=0;      % 2dim probability distributions
%close all % (de)activate for new/overwriting figures

% Data used [careful: in old version class was 'char' now it's 'double']
field1=dC5_cycCor;  % production rate 1. will be normalized to mean=1
field1str='dC5_cycCor'; % as string 
field2=dY5_cycCor;  % production rate 2. will be normalized to mean=1
field2str='dY5_cycCor'; % as string
field3=mu_Third_cycCor; % GROWTH RATE. all vectors need to have same length!
field3str='mu_Third_cycCor'; % as string
%***  FIELD 3 CAN BE NORMALIZED FACULATIVELY BELOW!!! ***

%field1=1+nG; % debugging
%field2=B;
%field3=1+mu;

N_gridpoints=2^8;             % # gridpoints for KDE (Power of 2. use 2^8 as default.)
Bandwidth_scaleField=1;     % sigma (bandwidth) of gaussian (for KDE): 
                            % automatic_sigma*Bandwidth_scaleField
                            % (automatic_sigma varies with the field)
                            % 'default' or 1: uses automaticaly optimized bandwidth
% Not implemented yet: allow different scale factors for different fields
cutoff_prob=0.001;          % cutoff below which conditional expectations are not plotted (too high errors)


% *********BLUBB FOR DEBUGGING ************
% make sure via (:)' that they are row vectors
%field1=1+rnd';
%field2=1+murnd2000(:)';%field1;%muRnd'.^2+yrndpart'.^2;
%field3=1+murnd2000(:)';

%;% murndsmall';
%field1str='YFPrnd';
%field2str='CFPrnd';
%field3str='Murnd';
%myschnitzname='corr =0.';
% % *********BLUBB END ************

% ************************************************************************
%% 2) Calculate explained (co)variances via KDE
% ************************************************************************

% make sure they're column vectors
field1=field1(:);
field2=field2(:);
%field3=field3(:);

% normalize production rates
field1=field1/mean(field1);
field2=field2/mean(field2);
field3=field3/mean(field3); %and growth rate?

% ****************************************
% ** Default output of called function: **
% [bandwidth_vecAvecB,Probjoint_vecAvecB, meshX_vecAvecB, meshY_vecAvecB, ...
%    Prob_vecA,Prob_vecB, Prob_vecA_given_vecB, ...
%    vecA_grid, vecB_grid, Mean_vecA_given_vecB, Varunexplained_vecA_knowing_vecB, ...
%    Varexplained_vecA_knowing_vecB, Vartotal_vecA_viaProbdistr,increment_dvecB]= ...
%    NW_Probabilities_Grid_viaKDE2(vecA,vecB,n_gridpoints,bandwidth_scaleA,bandwidth_scaleB)
% ****************************************
% field3 is the growth rate!

% ------------------------------------
% Extrinsic Noise
% ------------------------------------

% (A) Extrinsic noise: field1 & field3:
disp(' ')
disp(['Calculating ' field1str ' & ' field3str ':']);
 [bandwidth_f1f3,Probjoint_f1f3, meshX_f1f3, meshY_f1f3, Prob_f1,Prob_f3_fromf1, Prob_f1_given_f3, ...
    f1_grid, f3_grid_fromf1, Mean_f1_given_f3, ~, ~, ~,increment_df3_fromf1]= ...
    NW_Probabilities_Grid_viaKDE2(field1,field3,N_gridpoints,Bandwidth_scaleField,Bandwidth_scaleField);
bandwidth_f1=bandwidth_f1f3(1);
bandwidth_f3_fromf1=bandwidth_f1f3(2);

% (B) Extrinsic noise: field2 & field3:
disp(' ')
disp(['Calculating ' field2str ' & ' field3str ':']);
[bandwidth_f2f3,Probjoint_f2f3, meshX_f2f3, meshY_f2f3, Prob_f2,Prob_f3_fromf2, Prob_f2_given_f3, ...
    f2_grid, f3_grid_fromf2, Mean_f2_given_f3, ~, ~, ~,increment_df3_fromf2]= ...
    NW_Probabilities_Grid_viaKDE2(field2,field3,N_gridpoints,Bandwidth_scaleField,Bandwidth_scaleField);
bandwidth_f2=bandwidth_f2f3(1);
bandwidth_f3_fromf2=bandwidth_f2f3(2);

% (C) check if output for growth rates (field3) is consistent from both function calls
% (a) grid 
if min(f3_grid_fromf1-f3_grid_fromf2)~=0 | max(f3_grid_fromf1-f3_grid_fromf2)~=0 | increment_df3_fromf1~=increment_df3_fromf2
    error('Error: Grid spacing of field3 is not identical for both function calls!')
else
    f3_grid=f3_grid_fromf1;
    increment_df3=increment_df3_fromf1;
    clear f3_grid_fromf1 f3_grid_fromf2 increment_df3_fromf1 increment_df3_fromf2
end
% (b) probability distribution
if PLOTMUDISTRIBUTION
    mufig=figure;
    clf; hold on;
    set(gcf,'WindowStyle','docked')
    plot(f3_grid, Prob_f3_fromf1*increment_df3,'.b')
    plot(f3_grid, Prob_f3_fromf2*increment_df3,'.r')
    xlabel([field3str],'Interpreter','None')
    ylabel('Probability')
    legend(['field1: ' field1str],['field2: ' field2str],'location','NE')
    title([myschnitzname ': Different marginal prob distr of ' field3str ' for joint prob distr from KDE2D (fieldX,field3).'])
end
% (c) bandwidth
if bandwidth_f3_fromf1~=bandwidth_f3_fromf2
    disp(['WARNING: Automatic bandwidth of field3 (' field3str ') differs: '])
    disp(['(field1,field3): bw=' num2str(bandwidth_f3_fromf1) '.'])
    disp(['(field2,field3): bw=' num2str(bandwidth_f3_fromf2) '.'])
end
    
% for field3 (growth rate) use the averaged probability distribution as
% marginal probability distribution
Prob_f3=0.5*(Prob_f3_fromf1+Prob_f3_fromf2);
if PLOTMUDISTRIBUTION
    figure(mufig)
    plot(f3_grid, Prob_f3*increment_df3,'k','LineWidth',2)
    legend(['field1: ' field1str],['field2: ' field2str],'average','location','NE')
end

% (D) extrinsic noise: covariance decomposition
% total covariance Cov(f1,f2) (total extrinsic noise) -> uses the dataset,
% not infered probdistributions (what is better?)
Cov_f1f2=cov(field1,field2); % (N-1) normalization
Cov_f1f2=Cov_f1f2(1,2);

% explained covariance Cov(<f1|f3>,<f2,f3>) with prob distribution Prob(f3)
mean_f1=sum(Mean_f1_given_f3.*Prob_f3*increment_df3); % should be =1;
mean_f2=sum(Mean_f2_given_f3.*Prob_f3*increment_df3);
% <field1|field3>-<<field1>>
Relative_field1_given_field3=Mean_f1_given_f3-mean_f1;
Relative_field2_given_field3=Mean_f2_given_f3-mean_f2;
% explained covariance . which normalization? Well Sum(Probabilities)=1
Cov_explained_f1f2=sum(Relative_field1_given_field3.*...
    Relative_field2_given_field3.*Prob_f3*increment_df3);
% Fraction
FractionCov_explained_f1f2=Cov_explained_f1f2/Cov_f1f2;

% ------------------------------------
% Intrinsic Noise
% ------------------------------------
%
% Variance of (field1-field2)
diff_f1f2=field1-field2;
[bandwidth_diff,Probjoint_difff3, meshX_difff3, meshY_difff3, ...
   Prob_diff,Prob_f3_from_diff, ~, diff_grid, f3_grid_from_diff, Mean_diff_given_f3, ... 
   ~, Varexplained_diff_knowing_f3, ~,~]= ...
    NW_Probabilities_Grid_viaKDE2(diff_f1f2,field3,N_gridpoints,Bandwidth_scaleField,Bandwidth_scaleField);
% Total Variance
Var_f1minf2=var(field1-field2);
% Explained variance
Var_explained_f1minf2=Varexplained_diff_knowing_f3;
%Fraction
FractionVar_explained_f1minf2=Var_explained_f1minf2/Var_f1minf2;

% ------------------------------------------
% Report the results. The actual calculation is finished here
% ------------------------------------------

disp(' ')
disp('_______________________________________________________________')
disp(['Analyzed ' myschnitzname '.'])
disp('_______________________________________________________________')
disp(['Used KDE2D to determine probability distributions.'])
disp(['Used fields ' field1str ' & ' field2str ])
disp(['and conditioned on ' field3str ' without history.'])
disp('  ')
disp('EXTRINSIC')
disp('(BW=bandwidth)')
disp(['BW(mu)    factor*autoBW     #grid points   Explained_Extr      TotalNoise      Frac_Expl'])
disp([num2str(mean([bandwidth_f3_fromf1,bandwidth_f3_fromf2])) '           '  num2str(Bandwidth_scaleField) ...
    '           ' num2str(N_gridpoints) '        '  num2str(Cov_explained_f1f2) '              ' ...
    num2str(Cov_f1f2) '       '    num2str(FractionCov_explained_f1f2)])

disp('  ')
disp('INTRINSIC')
disp('(BW=bandwidth)')
disp(['BW(mu)    factor*autoBW     #grid points   Explained_Intr      TotalNoise      Frac_Expl'])
disp([num2str(bandwidth_diff(2)) '           '  num2str(Bandwidth_scaleField) ...
    '           ' num2str(N_gridpoints) '        '  num2str(Var_explained_f1minf2) '              ' ...
    num2str(Var_f1minf2) '       '    num2str(FractionVar_explained_f1minf2)])


% ************************************************************************
%% 3) Some Plotting
% ************************************************************************

if PLOTCONDVALUESEXTR
    idxplot=find(Prob_f3*increment_df3>cutoff_prob);
    
    figure
    clf; hold on;
    set(gcf,'WindowStyle','docked')
    subplot(2,1,1)
    hold on
    plot(f3_grid(idxplot), Mean_f1_given_f3(idxplot),'-','Color',[0 0.5 1],'LineWidth',2)
    plot(f3_grid(idxplot), Mean_f2_given_f3(idxplot),'-','Color',[1 0.5 0],'LineWidth',2)
    %title(['Conditional Expectations -> Extrinsic Noise. cutoff(Prob(mu)*dmu)='  num2str(cutoff_prob)]);
    title(['Cond Expect''s -> Extr Noise. cutoff(Prob(mu)*dmu)='  num2str(cutoff_prob)]);
    if (mean(field3)-1)<1e-10
        xlabel(['normalized ' field3str],'Interpreter','None')
    else
        xlabel([field3str],'Interpreter','None')
    end        
    ylabel('<Production|mu>')
    legend({['< ' field1str ' |mu>'],['< ' field2str ' |mu>']},'Interpreter','none','location','SE' )
    
    subplot(2,1,2)
    plot(f3_grid(idxplot), Prob_f3(idxplot),'-k','LineWidth',2)
    if (mean(field3)-1)<1e-10
        xlabel(['normalized ' field3str],'Interpreter','None')
    else
        xlabel([field3str],'Interpreter','None')
    end 
    ylabel('Probability(mu)')
    legend({'Prob(mu)'},'Interpreter','none','location','NE' )
    
end

if PLOTCONDVALUESINTR
    idxplot=find(Prob_f3_from_diff*increment_df3>cutoff_prob);
    
    figure
    clf; hold on;
    set(gcf,'WindowStyle','docked')
    plot(f3_grid(idxplot), Mean_diff_given_f3(idxplot),'-','Color',[0.3 0.3 0.3],'LineWidth',2)
   % plot(f3_grid(idxplot), 0.2*Prob_f3(idxplot)/max(Prob_f3(idxplot)),'.k','LineWidth',2)
    if (mean(field3)-1)<1e-10
        xlabel(['normalized ' field3str],'Interpreter','None')
    else
        xlabel([field3str],'Interpreter','None')
    end      
  %  ylabel('<field1-field2|mu> resp (rescaled) Probability')
    ylabel(['<' field1str '-' field2str '|mu>'])
    legend({['< ' field1str '-' field2str ' |mu>'],'Prob(mu)'},'Interpreter','none','location','NE' )
    %title(['Conditional Expectations -> Intrinsic Noise. cutoff(Prob(mu)*dmu)='  num2str(cutoff_prob)]);
    title(['Cond Expect''s -> Intr. Noise. cutoff(Prob(mu)*dmu)='  num2str(cutoff_prob)]);
end

%%
%if PLOTPROBHEATMAPS
%    BLUBB
%    
%    tt=Prob_f1_given_f3;
%    tt(Probjoint_f1f3*increment_df3<0.001)=NaN;
%surf(meshX_f1f3,meshY_f1f3,tt,'EdgeColor','None'); view([0 0 90])
%    
%    
%     surf(meshX_f1f3,meshY_f1f3,Probjoint_f1f3,'EdgeColor','None'); view([0 0 90])
%     
%     surf(meshX_difff3,meshY_difff3,Probjoint_difff3,'EdgeColor','None'); view([0 0 90])
    