function [bandwidth_vecAvecB,Probjoint_vecAvecB, meshX_vecAvecB, meshY_vecAvecB, ...
    Prob_vecA,Prob_vecB, Prob_vecA_given_vecB, ...
    vecA_grid, vecB_grid, Mean_vecA_given_vecB, Varunexplained_vecA_knowing_vecB, ...
    Varexplained_vecA_knowing_vecB, Vartotal_vecA_viaProbdistr, increment_dvecB]= ...
    NW_Probabilities_Grid_viaKDE2(vecA,vecB,n_gridpoints,bandwidth_scaleA,bandwidth_scaleB)
% This function calculates joint/marginal/conditional probabilities & Variances & 
% Expectation values for a bivariate dataset and also returns the coordinates of the 
% underlying grid (axis of 1st and 2nd vector).
% Input data is currently not normalized.
% There is more output in the direction vecA explained by vecB, "vecA
% conditioned on vecB" (e.g. CFP conditioned on mu)
%
% Its core is kde2d (2dimensional Kernel Density Estimate) (Matlab
% FileExchange ID #17204) which was slightly modified by adding extra
% options (leaving the essence untouched). kde2d estimates the joint
% probability distribution of a bivariate dataset by using a gaussian
% kernel.
% 
% The output of NW_Probabilities_Coordinates_viaKDE2.m is used to perform
% covariance (variance) decomposition for the extrinsic noise project.
%
% An alternative to calculate variance decomposition is ExtrinsicNoise_ExplaineByMu_VectorInput.m
% (less elaborate since binning of a finite dataset can cause
% overestimation of explained variance (probably neglegible for e.g. 15
% bins for 3000datapoints)
%
% NW 2014-12
%
% **************************
% INPUT
% **************************
% dimsensions are in square brackets
%
% vecA:             [n,1] 1st input vector (e.g. production rate CFP)
% vecB:             [n,1] 2nd input vector (e.g. growth rate mu)
% n_gridpoints:     size of the n by n grid over which the density is computed
%                   n has to be a power of 2, otherwise n=2^ceil(log2(n));
%                   If unsure, use 2^8; [directly fed intokde2d]
% bandwidth_scaleA: extra smoothing factor for bandwidth. Automatically
%                   calculated bandwidth for vecA is broadened by this
%                   factor (i.e. stddev of the gaussian kernel multiplied
%                   by this factor??)
%                   'default' or 1: use automatic bandwidth
%                   >1: extra smoothing
%                   <1: narrower gaussians
% bandwidth_scaleB:  same as bandwidth_scaleA for vecB
%
% cutoff          : This value is hardcoded below [i.e. no function input]! If probabilities in
%                   joint Probdistribution are <cutoff they are set to=0
%
% NOTE: kde2d allows to set the endpoints of the grid (MIN_XY, MAX_XY).
% This is not implemented as input argument in NW_Probabilities_Coordinates_viaKDE2
% and has to be set within this code if wished
%
% NOTE: This function provides some consistency check output (total
% probability etc). Activate by setting the hardcode TEXTOUTPUT=1
%
%
% **************************
% OUTPUT
% **************************
% dimensions are in square brackets
%
% bandwidth_vecAvecB:       [1,2] bandwidth for KDE: sigma auf gaussian (sigma_vecA,sigma_vecB)
% Probjoint_vecAvecB:       [n,n] joint probability distribution P(vecA,vecB)
% meshX_vecAvecB:           [n,n] meshgrid X coordinates for e.g. surf(meshX_vecAvecB, meshY_vecAvecB, Probjoint_vecAvecB)
% meshY_vecAvecB:           [n,n] meshgrid Y coordinates
% Prob_vecA,Prob_vecB:      [n,1] marginal probability distribution P(vecA), P(vecB)
% Prob_vecA_given_vecB:     [n,n] conditional probabilities P(vecA|vecB)
% vecA_grid:                [n,1] vecA coordinates for which prob distr is calculated
% vecB_grid:                [n,1] vecB ...
% Mean_vecA_given_vecB:     [n,1] <vecA|vecB>
% Varunexplained_vecA_knowing_vecB: [1] <Var(vecA|vecB)>
% Varexplained_vecA_knowing_vecB:   [1] Var(<vecA|vecB>)
% Vartotal_vecA_viaProbdistr:       [1] sum(probjoint*vecA^2*increment)-mean...
% increment_dvecB:          [1] sum(Prob_vecB)*increment_dvecB=1 (ouput for
%                                vecA currently not stored
%
% **************************
% Example Input
% **************************
% murnd=randn(2000,1);
% crnd2000_05=sqrt(1-0.5^2)*randn(2000,1)+0.5*murnd2000;
% NW_Probabilities_Grid_viaKDE2(crnd2000_05,murnd2000,2^8,'default','default')

TEXTOUTPUT=0;

% Descriptions of vectors and matrices are given wrt to (vecA,vecB) but also
% the typical inputs (CFP,mu)

% make sure vectors are column vectors
vecA=vecA(:);
vecB=vecB(:);
% normalize
%vecA=vecA/mean(vecA); 
%vecB=vecB/mean(vecB);

% check for default bandwidth
if strcmp(lower(bandwidth_scaleA),'default')==1
    bandwidth_scaleA=1;
end
if strcmp(lower(bandwidth_scaleB),'default')==1
    bandwidth_scaleB=1;
end
if ~isa(bandwidth_scaleA,'double') & ~isa(bandwidth_scaleA,'int')
    error('Error in input: bandwidth_scaleA must be scalar or the string ''default''.')
end
if ~isa(bandwidth_scaleB,'double') & ~isa(bandwidth_scaleB,'int')
    error('Error in input: bandwidth_scaleB must be scalar or the string ''default''.')
end

% ----------------
% KDE(2D)
% ----------------
% get joint probability distribution & coordinates of grid 
[bandwidth_vecAvecB,Probjoint_vecAvecB,meshX_vecAvecB,meshY_vecAvecB]= ...
        kde2d_inclbandwidth([vecA, vecB],n_gridpoints,bandwidth_scaleA,bandwidth_scaleB);
%kde2d([dC5_cycCor',mu_Third_cycCor'],2^8)

% example: kde2d([CFP, mu],...)
% X_vecAvecB, Y_vecAvecB: meshgrids

% ----------------
% grid coordinates
% ----------------
% get coordinates of probability distr grid: vectors with size (n_gridpoints,1)
vecA_grid=meshX_vecAvecB(1,:); % all A values for which prob. density is calculated (defult: 256) (e.g. CFP axis)
vecB_grid=meshY_vecAvecB(:,1); % all B values for which prob. density is calculated (defult: 256) (e.g. mu axis)
%old names: Caxis..., muaxis...
vecA_grid=vecA_grid(:);  % columns
vecB_grid=vecB_grid(:);

% ----------------
% increments dX
% ----------------
% increment of the axis -> for integration: Integral(P(A,B)*dA*dB)=1
% discrete integration: Sum(P(vecA,vecB)*increment_dvecA*increment_dvecB)=1
% (ideally sum=1)
increment_dvecA=mean(diff(vecA_grid)); % all diffs should be the same (equidistant)
increment_dvecB=mean(diff(vecB_grid));

% ----------------
% marginal probabilities
% ----------------
% Probjoint_vecAvecB is the joint probability P(vecB,vecA) ((mu,CFP)).
% each row of Probjoint_vecAvecB corresponds to a fixed vecB (mu) and runs through all
% vecA (CFP)
% i.e. row<->vecB and column<->vecA       ->Prob(row_vecB,column_vecA)
% 
%                      ( P(vecB1,vecA1)   P(vecB1,vecA2)   P(vecB1,vecA3),    ...   )
% Probjoint_vecAvecB=  ( P(vecB2,vecA1)   P(vecB2,vecA2)   P(vecB2,vecA3),    ...   )
%                      ( P(vecB3,vecA1)   ...             ...
%                           ...
%
%
% For the CFP&mu case: row<->mu and column<->CFP       ->Prob(row_mu,column_CFP)
% 
%                   ( P(mu1,CFP1)   P(mu1,CFP2)   P(mu1,CFP3),    ...   )
% Probjoint_CFPmu=  ( P(mu2,CFP1)   P(mu2,CFP2)   P(mu2,CFP3),    ...   )
%                   ( P(mu3,CFP1)   ...             ...
%                       ...

% marginal distributions: Prob_vecA=Sum(P(vecA,vecB)*increment_dvecB (integrate over vecB)
% P(vecA): one needs to sum over all rows  (P(CFP))
% P(vecB): one needs to sum over all columns  (P(mu))
Prob_vecA=sum(Probjoint_vecAvecB,1)*increment_dvecB; % 1: sum along 1st dimension (rows)
Prob_vecB=sum(Probjoint_vecAvecB,2)*increment_dvecA;  % 2: sum along 2nd dimension (columns)
% convert them into/make sure that they are column vectors 
Prob_vecA=Prob_vecA(:);
Prob_vecB=Prob_vecB(:);

% ----------------
% cutoff probabilities
% ----------------
% The matrix Probjoint has some very low values (e.g. <10^-5). When
% calculating the conditional probabilities, divison by these low values
% leads to absurd spikes in the cond. prob. (divergence)
% 1) if wishing to plot e.g. conditional probabilities: define cutoff
% (depends on scale of input data! e.g. 0.0001)
% 2) for law of total (co)variance the spikes seem to cancel when outer
% integral (over marginal prob of mu) is performed) and best results are
% obtained when no cutoff is introduced
%
% alternatively the cutoff could be performed (worse results?) before
% determining marginal distributions

cutoff=-1;%0.001; %-1: nothing excluded
idx=find(Probjoint_vecAvecB<cutoff); 
Probjoint_vecAvecB(idx)=0;


% ----------------
% conditional probabilities
% ----------------
% determine the conditional distributions P(vecA|vecB) and P(vecB|vecA)   (P(CFP|mu) and P(mu|CFP))
% P(vecA|vecA)=P(vecB,vecA)/P(vecB)

% Prob_vecA_given_vecB: The first row of Probjoint_vecAvecB has to be
%                   divided by Prob(vecB=vecB1) (=the first entry of
%                   Prob_vecB
% in CFP-mu terms: Prob_CFP_given_mu: The first row of Probjoint_Cmu has to be 
%                   divided by Prob(mu1) (=the first entry of Prob_mu)
%
%
%                        ( P(vecA1 | vecB1)   P(vecA2 | vecB1)   P(vecA3 | vecB1),    ...   )
% Prob_vecA_given_vecB=  ( P(vecA1 | vecB2)   P(vecA2 | vecB2)   P(vecA3 | vecB2),    ...   )
%                        ( P(vecA1 | vecB3)   ...             ...
%                            ...
%
% in CFP-mu terms:
%                     ( P(CFP1 | mu1)   P(CFP2 | mu1)   P(CFP3 | mu1),    ...   )
% Prob_CFP_given_mu=  ( P(CFP1 | mu2)   P(CFP2 | mu2)   P(CFP3 | mu2),    ...   )
%                     ( P(CFP1 | mu3)   ...             ...
%                        ...

% elementwise division
Prob_vecA_given_vecB=bsxfun(@rdivide,Probjoint_vecAvecB,Prob_vecB); % divide by column vec
Prob_vecB_given_vecA=bsxfun(@rdivide,Probjoint_vecAvecB,Prob_vecA'); %divide by row vec
%Prob_CFP_given_mu=bsxfun(@rdivide,Probjoint_vecAvecB,Prob_mu_fromvecAvecB); % divide by column vec
%Prob_mu_given_CFP=bsxfun(@rdivide,Probjoint_vecAvecB,Prob_CFP_fromvecAvecB'); %divide by row vec

% if there are NaN entries (zero probability P(x,y)/P(y)=0/0 assign them to 0)
idxnan=find(isnan(Prob_vecA_given_vecB));
Prob_vecA_given_vecB(idxnan)=0;
idxnan=find(isnan(Prob_vecB_given_vecA));
Prob_vecB_given_vecA(idxnan)=0;

% ** debugging start **
% checking the sum (can only be correct if items not set to=0)
% <1 happens for areas of low total probabilty (=at the edges). but these
% areas are also not important for calculating total expectation values etc
% figure
% hist(sum(Prob_vecB_given_vecA,1)*increment_dvecB,100)
% title('integrated conditional probabilities')
% xlabel('sum(Prob_vecB_given_vecA)*increment_dvecB -> should be=1 (chance
% of any B = 1)')
% ylabel('frequency')
% figure
% hist(sum(Prob_vecA_given_vecB,1)*increment_dvecA,100)
% title('integrated conditional probabilities')
% xlabel('sum(Prob_vecA_given_vecB)*increment_dvecA -> should be=1 (chance
% of any A = 1)')
% ylabel('frequency')

%figure(1)
%surf(X_vecAvecB,Y_vecAvecB, Prob_vecB_given_vecA,'EdgeColor','None')
%xlabel('vecA');ylabel('vecB')

% ** debugging end **

% 
% ----------------
% chech conservation of probabilities
% ----------------
% total probability seems to be >1. maybe include a correction factor?
totalProbability=sum(sum(Probjoint_vecAvecB))*increment_dvecA*increment_dvecB;
minProb=min(min(Probjoint_vecAvecB));
maxProb=max(max(Probjoint_vecAvecB));
totalProb_A=sum(Prob_vecA)*increment_dvecA;
totalProb_B=sum(Prob_vecB)*increment_dvecB;
if TEXTOUTPUT
    disp(' ')
    disp('..........')
    disp('Total probabilities (check for conservation of probability -> sum=1) :')
    disp(['sum(joint probabilities P(vecB,vecA)) = ' num2str(totalProbability)])
    disp(['minimum joint prob. density (<0?, magnitude compared to max?): ' num2str(minProb)])
    disp(['max joint prob. density: ' num2str(maxProb) '. '])
    disp(['Total marginal probability sum(Prob(vecA)): ' num2str(totalProb_A)]);
    disp(['Total marginal probability sum(Prob(vecB)): ' num2str(totalProb_B)]);
    disp('..........')
end


% ----------------
% conditional expectation values
% ----------------
% <vecA|vecB1>: sum the first row of density matrix * vecA_grid * increment_dvecA
% in CFP&mu terms: <CFP|mu1>  : sum the first row of density matrix * Caxis * increment_dC
Mean_vecA_given_vecB_temp=bsxfun(@times,Prob_vecA_given_vecB,vecA_grid(:)')*increment_dvecA; % row vector!
Mean_vecA_given_vecB=sum(Mean_vecA_given_vecB_temp,2); % each entry: a different vecB_i
%MeanCFP_given_mu_temp=bsxfun(@times,Prob_CFP_given_mu,Caxis_fromCmu(:)')*increment_dC; % row vector!
%MeanCFP_given_mu=sum(MeanCFP_given_mu_temp,2);
Mean_vecA_given_vecB=Mean_vecA_given_vecB(:);

% % inverted conditional expectation: to be tested!
% % <vecB|vecA1>
%Mean_vecB_given_vecA_temp=bsxfun(@times,Prob_vecB_given_vecA,vecB_grid(:))*increment_dvecB; % column vector!
%Mean_vecB_given_vecA=sum(Mean_vecB_given_vecA_temp,1); % each entry: a different vecA_i
%Mean_vecB_given_vecA=Mean_vecB_given_vecA(:);


% ----------------
% (mean) conditional variance. blubb to be tested
% ----------------
% for intrinsic noise (unexplained variance):
% <Var(vecA|vecB)>   [ inverted conditioning (vecB|vecA) is not calculated]
% in CFP/mu terms: <Var(CFP|mu> or <Var(CFP-YFP)|mu>
Var_vecA_given_vecB_temp=bsxfun(@times,Prob_vecA_given_vecB,vecA_grid(:)'.*vecA_grid(:)')*increment_dvecA; % column vector!
Var_vecA_given_vecB=sum(Var_vecA_given_vecB_temp,2);
Var_vecA_given_vecB=Var_vecA_given_vecB-Mean_vecA_given_vecB.^2; % column vec: [Var(vecA|B1); Var(vecA|B2); ...]
Varunexplained_vecA_knowing_vecB=sum(Var_vecA_given_vecB.*Prob_vecB*increment_dvecB); % scalar

% ----------------
% variance of conditional expectations. blubb to be tested
% ----------------
% for intrinsic noise (explained variance)
% Var(<vecA|vecB>)
Varexplained_vecA_knowing_vecB=sum(increment_dvecB*Prob_vecB.*Mean_vecA_given_vecB.*Mean_vecA_given_vecB) ...
    -(sum(increment_dvecB*Prob_vecB.*Mean_vecA_given_vecB))^2;


% ---------------
% total variance Var(vecA)
% ---------------
% Var(vecA)
% from the input vector directly
Vartotal_vecA_viaInputVector=var(vecA);
% via calculated Prob distribution (should be the same!)
Vartotal_vecA_viaProbdistr=sum(Prob_vecA.*vecA_grid.^2*increment_dvecA)-...
    sum(Prob_vecA.*vecA_grid*increment_dvecA)^2;

%  UNTIL HERE KIND OF DONE
 % explained variance (to do: covariance blubb)
%Varxepl=sum(increment_dvecB*Prob_vecB.*Mean_vecA_given_vecB.*Mean_vecA_given_vecB) ...
%    -(sum(increment_dvecB*Prob_vecB.*Mean_vecA_given_vecB))^2

% %Vartotal=
%disp(['first: ' num2str(sum(increment_dvecB*Prob_vecB.*Mean_vecA_given_vecB.*Mean_vecA_given_vecB)) ...
%    '   2nd:  ' num2str((sum(increment_dvecB*Prob_vecB.*Mean_vecA_given_vecB))^2)...
%    '      mean (sqrt of 2nd):'
%    num2str(sum(increment_dvecB*Prob_vecB.*Mean_vecA_given_vecB))]);%

if TEXTOUTPUT
    disp(['varexpl: ' num2str(Varexplained_vecA_knowing_vecB) ...
        '   varnotexpl:  ' num2str(Varunexplained_vecA_knowing_vecB)...
        '   vartotal_viaprobdistr:  ' num2str(Vartotal_vecA_viaProbdistr)...
        '      vartotal_viainput: ' num2str(Vartotal_vecA_viaInputVector)]);
    disp(['frac explained (div. by Var_viainputvector):  ' num2str(Varexplained_vecA_knowing_vecB/Vartotal_vecA_viaInputVector)])
    disp(' ')
end

% NB: Total variance of the Probdistribution differs somewhat from total
% variance of the input vector!
% For Cov, the probdistr would have to be calculated of yet another pair
% (vecA,vecB)=(YFP,CFP) -> maybe simply use the experimental distr?


% do the axis vectors from different kde2d have hte same values? should be
% because i think they are ecided solely on the inputdatavector

% potential problem: the marginal distributions P(mu) which are obtained b
% "integrating" P(CFP,mu) over CFP resp P(YFP,mu) over YFP are somewhat
% different. which one to choose? or calculate it via a univariate KDE?
