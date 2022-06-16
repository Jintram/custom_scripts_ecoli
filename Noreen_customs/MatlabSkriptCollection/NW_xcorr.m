function [mycrosscorr, timelags] = NW_xcorr(datavec1,datavec2,varargin)
% Computes the cross correlation (respectively autocorrelation if input
% vectors are equal)
% The function is based on the matlab function 'xcorr'. However, several
% normalizing features are implemented as default setting or added.
%
% -------------------------------------------------------------
% Default setting: Cross correlation is computed according to:
% -------------------------------------------------------------
%
% CC(tau)=[<x(t)*y(t+tau)> - <x>*<y>] / [ sttdev(x) * stddev(y) ]
%
% The matlab funxtion xcorr displays as default result: Sum[x(t)*y(t+tau)].
% Here, several features are added as default: 
% 1) division by the number of data points ('unbiased' in xcorr) to obtain the actual 
%    average square <x(t)*y(t+tau)>   ( this 'unbiased' normalization is
%    important because at large time delays, the amount of data is less and
%    this has to be corrected for (i.e. divide by roughly (N-timelag)
%    (N=length of data vector)
% 2) subtraction of the mean <x><y> since we are mostly interested in
%    fluctuations around a mean
% 3) Normalization by the standard deviations stddev(x)*stddev(y). This
%    transforms the crosscovariance into a crosscorrelation (meaning that
%    the autocorrelation is normalized to =1 at time lag 0)
%
% Options 2) and 3) can be altered via extra inputs
%         
% -------------------------------------------------------------
% Input
% -------------------------------------------------------------
% datavec1: first input vector x for crosscorrelation (length N)
% datavec2: second input vector y for crosscorrelation (length N). If datavec1 and
%           datavec2 are identical, the autocorrelation will be computed
% **varargin**
% 'timeincr': optional. increment of the time datapoints. preceded by the argument
%            'timeincr'  (...,'timevec',somedouble,...)
% 'meansubtr': optional. default=1 -> <x><y> is subtracted. If set to =0,
%           mean is not subtracted (...'meansubtr',0,...)
% 'stddevnorm' : optional. default=1 -> normalization with stddev(x)*stddev(y).
%                if set to =0, no extra normalization is performed.
%                (...,'stddevnorm',0...)
%
% -------------------------------------------------------------
% Output
% -------------------------------------------------------------
% mycrosscorr: cross correlation vector of size (2N-1). If autocorrelation:
%              size is (N)
% timelags:    vector of size (2N-1) containing all corresponding time
%              lags. If no time vector was given as input the timelags are
%              default steps [-N+1:1:N-1]. If autocorrelation: vector of size
%              (N) and in abscence of a time vector: [0:1:N-1]
%
% -------------------------------------------------------------
% example function call: [cc,time] = NW_xcorr(X,Y,'timeincr',0.001,'meansubtr',0)
% ------------------------------------------------------------------------------------------


% -------------------------------------------------------------
% Set Default Values and Get Function Inputs (varargin) (overwrite default)
% -------------------------------------------------------------
timeincr=1;
meansubtr=1;
stddevnorm=1;

if size(varargin,1)>0 % extra arguments given
    if mod(size(varargin,2),2)~=0
        error('Extra arguments must be given in pairs.')
    else
        for i=1:size(varargin,2)-1
            argum=varargin{i};
            if strcmp(argum,'timeincr')==1
                timeincr=varargin{i+1};
            elseif strcmp(argum,'meansubtr')==1
                meansubtr=varargin{i+1};
            elseif strcmp(argum,'stddevnorm')==1
                stddevnorm=varargin{i+1};
            end
        end
    end
end
                
% check for some input errors
if length(datavec1)~=length(datavec2)
    error('Input vectors do not have the same length')
end
if meansubtr~=0 & meansubtr~=1
    disp('Wrong ''meansubtr'' input. Will use =1.')
    meansubtr=1;
end
if stddevnorm~=0 & stddevnorm~=1
    disp('Wrong ''stddevnorm'' input. Will use =1.')
    stddevnorm=1;
end

% prepare time vector
timevec=[1:length(datavec1)]*timeincr;
    
% -----------------------------------------------------------------------------
% Compute cross correlation
% -----------------------------------------------------------------------------
xy_mean=xcorr(datavec1,datavec2,'unbiased');
mycrosscorr=xy_mean;

if meansubtr==1
    mycrosscorr=mycrosscorr-mean(datavec1)*mean(datavec2);
end
if stddevnorm==1
    mycrosscorr=mycrosscorr/(std(datavec1,1)*std(datavec2,1)); % use flag=1 -> normalize by N and not (N-1)
end

timelags=[-timevec(end-1:-1:1),0,timevec(1:end-1)];

% check if input vectors are identical -> autocorrelation -> restrict to
% times >=0
difference1_2=sum(abs(datavec1-datavec2));
if difference1_2==0
    disp('Input vectors are identical. Will compute autocorrelation.')
    N=round(length(mycrosscorr)-1)/2; % size of initial input vectors
    mycrosscorr=mycrosscorr(end-N:end);
    timelags=timelags(end-N:end);
end

