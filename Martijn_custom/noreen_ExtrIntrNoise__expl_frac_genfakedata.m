%%

N = 3512
FRACTION = .2
gainfactor = 2
Offset = 0; % just to get positive mu-values

% normal distributions, mu = growth rate; Y and C are fluor signals
muRnd       = randn(N,1)+Offset;
yrndpart    = randn(N,1)+Offset; % intrinsic noise in Y
crndpart    = randn(N,1)+Offset; % intrinsic noise in C

% 
YFPrnd = gainfactor.* (FRACTION .* muRnd + sqrt(1-FRACTION^2) .* yrndpart);
CFPrnd =               FRACTION .* muRnd + sqrt(1-FRACTION^2) .* crndpart;

%YFPrnd=0.2535*muRnd+sqrt(1-0.2535^2)*yrndpart;
%CFPrnd=0.1520*muRnd+0.6*yrndpart+(1-0.6^2-0.152^2)*crndpart;

figure(1); clf;
plot(muRnd, YFPrnd, '.r')
hold on;
plot(muRnd, CFPrnd, '.b')

% Note that there is no time component in this data, so we cannot make
% cross correlation functions.
% Only the correlation coefficient can be measured.

RmuY = corr(muRnd,YFPrnd)
RmuC = corr(muRnd,CFPrnd)

% Note that Noreen has sent me two websites with two nice explanations when
% the slope indeed is the same as the cov, and when not.
%
% - http://stats.stackexchange.com/questions/22718/what-is-the-difference-between-linear-regression-on-y-with-x-and-x-with-y/22721#22721
% - http://stats.stackexchange.com/questions/32464/how-does-the-correlation-coefficient-differ-from-regression-slope

%%

%{

muRnd=randn(N,1);
yrndpart=randn(N,1);
crndpart=randn(N,1);
%YFPrnd=0.2535*muRnd+sqrt(1-0.2535^2)*yrndpart;
%CFPrnd=0.1520*muRnd+0.6*yrndpart+(1-0.6^2-0.152^2)*crndpart;
YFPrnd02=FRACTION*muRnd+sqrt(1-FRACTION^2)*yrndpart;
CFPrnd02=FRACTION*muRnd+sqrt(1-FRACTION^2)*crndpart;

muRnd=randn(3512,1);
yrndpart=randn(3512,1);
crndpart=randn(3512,1);
%YFPrnd=0.2535*muRnd+sqrt(1-0.2535^2)*yrndpart;
%CFPrnd=0.1520*muRnd+0.6*yrndpart+(1-0.6^2-0.152^2)*crndpart;
YFPrnd02=0.2*muRnd+sqrt(1-0.2^2)*yrndpart;
CFPrnd02=0.2*muRnd+sqrt(1-0.2^2)*crndpart;

%}

