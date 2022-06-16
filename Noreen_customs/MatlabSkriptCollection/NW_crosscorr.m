function [vecdelaytime, crosscorr, taudecay, exponentialfit]=NW_crosscorr(veca,vecb,varargin);
% This function calculates the cross-correlation of two input vectors veca
% & vecb. The time vector can be given as optional 3rd input. 
% Cross-corrs are mean-subtracted, normalized by stddev and unbiased (=take care of
% less data at larger delays). 
% If veca==vecb, i.e. computation of auto-correlation, the output vectors are truncated 
% to positive delays and an exponential decay function is fitted.

% ------------
% input
% ------------
% veca, vecb:     input vectors for cross-corr. same length. can be identical
% optinal (varargin): timevector

% ------------
% output
% ------------

% ----------------------------------------------------
% parse inputs. some error checking
% ----------------------------------------------------
% vecdelaytime:     [-max_timedlay ... 0 ... +max_timedlay], truncated to
%                   [0 ... + max_timedelay] for autocorrs
% crosscorr:        obvious. truncated to positive delays (from 0 on) for autocorrs
% taudecay:         ***exponential*** decay time (when is autocorr
%                    decreased to 1/e). only calculated for autocorr
% exponentialfit:   object with fit parameters. only calc for autocorr
% do some trivial input checks
if length(veca)~=length(vecb)
    error('Input vectors must be same length');
end

%make sure we work with column vectors
veca=veca(:);
vecb=vecb(:);

% check if time vector is given, otherwise use [1...N]'
if nargin==2;
    vectime=[1:1:length(veca)]';
else
    vectime=varargin{1};
    vectime=vectime(:);
    if ~size(vectime)==size(veca)
        error('Time vector has wrong length.');
    end
end
% check for equidistant time intervals
deltatime=mean(diff(vectime));
timedeviation=max(diff(vectime)-deltatime)/deltatime;
tolerancedeviation=0.0001; % threshold for double numbers for which time points are considered equidistant
if abs(timedeviation)>tolerancedeviation
    error('Time points not equidistant');
end

% #datapoints
numdata=length(veca);


% ----------------------------------------------------
% calculate crosscorr
% ----------------------------------------------------

% todo: get the +- sign right!!
% CC(a,b,tau)=[<a(t),b(t+tau)>-<a><b>] / [std(a)*std(b)]

% subtract mean
veca=veca-mean(veca);
vecb=vecb-mean(vecb);

% get stddev
std_veca=std(veca);
std_vecb=std(vecb);

% crosscorr. use inbuilt fct xcorr with 'unbiased' (scales the raw
%       correlation by 1/(M-abs(lags)) = takes care of less data at higher
%       delays)
crosscorr=xcorr(veca,vecb,'unbiased')/(std_veca*std_vecb);

% delay time vector
starttime=vectime(1);
postimedelays=vectime-starttime;
postimedelays=postimedelays(2:end); % exclude the delay=0 for the moment
negtimedelays=-postimedelays(end:-1:1);

vecdelaytime=[negtimedelays; 0; postimedelays];


% ----------------------------------------------------
% if autocorr: retrict to positive times and calculate decay time constant
% ----------------------------------------------------
taudecay=nan;
exponentialfit=nan;

% identical input vectors -> autocorr
if veca==vecb
    vecdelaytime=vecdelaytime(numdata:end);
    crosscorr=crosscorr(numdata:end);

    % fit the decay with an exponential function
    %       based on automatic code generation of matlab for an exponential fit (in
    %       the curve fitting session)
    %  Output:
    %      fitresult : a fit object representing the fit.
    %      gof : structure with goodness-of fit info.
    %
    [xData, yData] = prepareCurveData( vecdelaytime, crosscorr ); % check for dimensions, NaN etc
    % Set up fittype and options.
    ft = fittype( 'exp1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.0665820677876644 -0.28295353911268];   % apparently matlab default
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    % % Plot fit with data.
    % figure
    % h = plot( fitresult, xData, yData );

    exponentialfit=fitresult; %renaming the object
    taudecay=-1/exponentialfit.b;
end