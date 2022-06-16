% determine the explained fraciton with the help of KDE (kernel density
% estimate) instead of direct binning

% Therefore we need (at least) the probabilities
% P(YFP),P(CFP),P(mu),P(YFP,mu),P(CFP,mu)

% potential problem: the marginal distributions P(mu) which are obtained b
% "integrating" P(CFP,mu) over CFP resp P(YFP,mu) over YFP are somewhat
% different. which one to choose? or calculat it via a univariate KDE?

% This script is uses kde2d.m (from Matlab FileExchange ID #17204)

[bandwidth_Cmu,density_Cmu,X_Cmu,Y_Cmu]=kde2d([crnd2000_sqr, murnd2000],2^8);%kde2d([Crnd07, muRnd],2^8);%kde2d([CRnd,muRnd2]);%
%kde2d([dC5_cycCor',mu_Third_cycCor'],2^8)
% in the following CFP is first entry, mu is second entry!!!

Caxis_fromCmu=X_Cmu(1,:); % all cfp values for which prob. density is calculated (defult: 256)
                          % X_Cmu is a meshgrid
muaxis_fromCmu=Y_Cmu(:,1); % all mu values for which prob. density is calculated (defult: 256)

% increment of the axis -> for integration: Integral(P(y,z)*dY*dZ)=1
% (ideally)
% discrete integration: Sum(P(x,y)*increment_dC*increment_dmu)=1
increment_dC=mean(diff(Caxis_fromCmu)); % all diffs should be the same (equidistant)
increment_dmu=mean(diff(muaxis_fromCmu));

% density_Cmu is the joint probability P(mu,CFP)
% each row of density_Cmu corresponds to a fixed mu and runs through all
% CFP
% i.e. row<->mu and column<->CFP       ->Prob(row_mu,column_CFP)
% 
%               ( P(mu1,CFP1)   P(mu1,CFP2)   P(mu1,CFP3),    ...   )
% density_Cmu=  ( P(mu2,CFP1)   P(mu2,CFP2)   P(mu2,CFP3),    ...   )
%               ( P(mu3,CFP1)   ...             ...
%                   ...

% determining the marginal distributions
% P(mu): one needs to sum over all columns
% P(CFP): one needs to sum over all rows
Prob_mu_fromCmu=sum(density_Cmu,2)*increment_dC;  % 2: sum along 2nd dimension (columns)
Prob_CFP_fromCmu=sum(density_Cmu,1)*increment_dmu; % 1: sum along 1st dimension (rows)
% convert them into/make sure that they are column vectors 
Prob_CFP_fromCmu=Prob_CFP_fromCmu(:);
Prob_mu_fromCmu=Prob_mu_fromCmu(:);

% set very low probabilities to =0 (high errors! especially divergence in
% division!)  (do this after calculating/summing the marginal distributions


cutoff=-1;%0.001; %blubb adjust %for CFP & mu & 3500datapts: 0.001 was good. but scales with variables!!!! BLUBB
idx=find(density_Cmu<cutoff); % IMPORTANT!!!??
density_Cmu(idx)=0;


% determine the conditional distributions P(CFP|mu) and P(mu|CFP) (the
% latter is not needed but for the sake of completeness)
% P(CFP|mu)=P(mu,CFP)/P(mu)

% Prob_CFP_given_mu: The first row of density_Cmu has to be divided by
% Prob(mu1) (=the first entry of Prob_mu_fromCmu)

%                     ( P(CFP1 | mu1)   P(CFP2 | mu1)   P(CFP3 | mu1),    ...   )
% Prob_CFP_given_mu=  ( P(CFP1 | mu2)   P(CFP2 | mu2)   P(CFP3 | mu2),    ...   )
%                     ( P(CFP1 | mu3)   ...             ...
%                   ...

%Prob_CFP_fromCmu=density(:); % delete

% elementwise division
Prob_CFP_given_mu=bsxfun(@rdivide,density_Cmu,Prob_mu_fromCmu); % divide by column vec
Prob_mu_given_CFP=bsxfun(@rdivide,density_Cmu,Prob_CFP_fromCmu'); %divide by row vec

% checking the sum (can only be correct if items not set to=0)
%hist(sum(Prob_mu_given_CFP,1)*increment_dmu,100)
%hist(sum(Prob_CFP_given_mu,2)*increment_dC,100)

%figure(1)
%surf(X_Cmu,Y_Cmu, Prob_mu_given_CFP,'EdgeColor','None')
%xlabel('CFP');ylabel('mu')

% if there are NaN entries (zero probability P(x,y)/P(y)=0/0 assign them to
% 0)
idxnan=find(isnan(Prob_mu_given_CFP));
Prob_mu_given_CFP(idxnan)=0;
idxnan=find(isnan(Prob_CFP_given_mu));
Prob_CFP_given_mu(idxnan)=0;

%CHECK PROBABILITIES summing


% COND EXPCTS
% <CFP|mu1>  : sum the first row of density matrix * Caxis * increment_dC
MeanCFP_given_mu_temp=bsxfun(@times,Prob_CFP_given_mu,Caxis_fromCmu(:)')*increment_dC; % row vector!
%MeanCFP_given_mu_temp=bsxfun(@times,Prob_CFP_given_mu,ones(256,1))*increment_dC; %blubb

MeanCFP_given_mu=sum(MeanCFP_given_mu_temp,2);


% explained variance (to do: covariance blubb)
Varxepl=sum(increment_dmu*Prob_mu_fromCmu.*MeanCFP_given_mu.*MeanCFP_given_mu) ...
    -(sum(increment_dmu*Prob_mu_fromCmu.*MeanCFP_given_mu))^2

%Vartotal=
disp(['first: ' num2str(sum(increment_dmu*Prob_mu_fromCmu.*MeanCFP_given_mu.*MeanCFP_given_mu)) ...
    '   2nd:  ' num2str((sum(increment_dmu*Prob_mu_fromCmu.*MeanCFP_given_mu))^2)...
    '      mean (sqrt of 2nd):' num2str(sum(increment_dmu*Prob_mu_fromCmu.*MeanCFP_given_mu))]);



% do the axis vectors from different kde2d have hte same values? should be
% because i think they are ecided solely on the inputdatavector


%%
% sum over marginal probabilities * C_value
sum(Prob_CFP_fromCmu.*Caxis_fromCmu')*increment_dC

%should be the same as
sum(increment_dmu*Prob_mu_fromCmu.*MeanCFP_given_mu)
