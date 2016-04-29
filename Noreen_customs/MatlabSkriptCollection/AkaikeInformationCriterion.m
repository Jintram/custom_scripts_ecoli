function AIC = AkaikeInformationCriterion(sigma2_est,n,K,varargin)
% Implementation of the AKAIKE INFORMATION CRITERION (AIC), for further info,
% consult the book Burnham,Anderson: Model Selection and Multimodel
% Interference (2002).
% PROBLEM SETTING: 
% A dataset can be fitted with different models (e.g. a
% linear and a bilinear model). Often, the sum of the squared residuals
% (Mean square deviation) is used as criterion for goodness of fit.
% However, if one model contains more fitting parameters than another, the
% fit will most likely be better, up to overfitting. With very many fitting
% parameters a perfect fit can be produced.
% Thus, a penalty for the number of fitting parameters has to be introduced
% when judging the goodness of fit.
% A modern option to account for fitting parameters is the AIC
%
% *********************************************************************
% AIC = -2*log (L(Theta)) + 2*K + 2K(K+1)/(n-K-1)         (natural log)
%
% n: sample size
% K: number of estimated regression parameters (that is the free fitting
%    parameters plus(!) the variance), e.g.: linear fit y=ax+b  -> K=3
% L(Theta): Likelihood of the maximum-likelihood-estimate Theta.
%
% ---------- Lower AIC -> better model --------------
%
% *********************************************************************
% This function currently only works for NORMALLY DISTRIBTUED INDEPENDENT 
% ERRORS (fit deviations) that have a constant variance [ this probably
% implies a LINEAR LEAST SQUARE fit, I'm not sure about these details ]
% Then, the likelihood can be expressed via the root mean square deviation
% of the fitted data points:
% sigma^2 = 1/n * (Sumalldatapoints over (deviations^2))
% log(L(Theta))=-1/2*n*log(sigma^2)
%
% Then (formula used in this function): 
% -------------------------------------------------
% AIC = n*log (sigma^2) + 2*K + 2K(K+1)/(n-K-1)
% -------------------------------------------------
%
% Only relative AICs are relevant, therefore it is insensitive to rescaling
% of the dataset -> This would only introduce an extra constant in the AIC
% due to the log(sigma^2)
%
% *************************************************************************
%
% INPUT ARGUMENTS
% sigma2_est: estimated mean squared(!) deviation (error) from least square
%             fitting: sigma2_est=1/n Sum(i=1->n) (epsilon_i)^2
%             (epsilon=distance to fit = y-yfit). Don't forget the
%             normalization by n!
%             Other naming: MSE
% n: number of data points
% K: number of regression parameters: fitting parameters + variance (don't
%             forget the plus 1!)
% varargin: 'largesample' Then the last term of the AIC, which corrects for
%           small sample size, is omitted (no real point in using it)
%
% OUTPUT
% AIC: A value that estimates the quality of a model. Small value
%      -> good model
%      Only the relative differences of AICs of different models are
%      relevant. Constants are irrelevant if they appear for all models.


% set default AIC computation
LARGESAMPLE=0;

% check input
if ~isempty(varargin)
    checklargesample=varargin(1);
    if strcmp(checklargesample,'largesample')==1
        LARGESAMPLE=1;
    else
        disp('Input of varargin unknown. See help function.')
    end
end

if ~LARGESAMPLE
    AIC=n*log(sigma2_est) + 2*K + 2*K*(K+1)/(n-K-1);
else
    AIC=n*log(sigma2_est) + 2*K;
end
