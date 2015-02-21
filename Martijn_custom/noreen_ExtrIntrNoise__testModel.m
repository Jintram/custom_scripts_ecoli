%%
% Written by MW, 2015/02
% Script test model Noreen drafted to explain intr/extr noise w. different
% noise sources working different observables, some observables sharing
% noise sources.

N = 999999

% There are 5 noise sources
amplitudeMu     = 2; % editable parameters
amplitudeY      = 1;
amplitudeC      = 3;
amplitudeYC     = 4;
amplitudeMuYC   = 5;

% TODO: offset?
noiseMu     = randn(N,1) * amplitudeMu; %
noiseY      = randn(N,1) * amplitudeY; % also known as intrinsic noise
noiseC      = randn(N,1) * amplitudeC; % idem
noiseYC     = randn(N,1) * amplitudeYC;
noiseMuYC   = randn(N,1) * amplitudeMuYC; % also known as extrinsic noise
    % TODO: actually, is noiseYC considered extr./intr.? and noiseMu?

% These have effects on three observables
% Effects are defined through gains
noiseMuGain      = 2.5; % editable parameters
noiseYGain       = 3;
noiseCGain       = 2;
noiseYCGainY     = 1.5;
noiseYCGainC     = 4;
noiseMuYCGainMu  = 5;
noiseMuYCGainY   = 4.5;
noiseMuYCGainC   = 3.5;

% Now calculate the signal for the observables:
observableY     = noiseYGain*noiseY     + noiseYCGainY*noiseYC + noiseMuYCGainY*noiseMuYC;
observableC     = noiseCGain*noiseC     + noiseYCGainC*noiseYC + noiseMuYCGainC*noiseMuYC;
observableMu    = noiseMuGain*noiseMu   +                        noiseMuYCGainMu*noiseMuYC;
    % TODO: NoiseMuY and noiseMuC don't exist because the promoters are
    % symmetrical. At least, that's an assumption -- which might be 
    % incorrect if you consider the data.

figure(1); clf;
plot(observableMu, observableY, '.r')
hold on;
plot(observableMu, observableC, '.b')
axis equal
xlabel('Growth speed \mu')
ylabel('Fluor signals (Y & C)')

figure(2); clf;
plot(observableY, observableC, '.r')
axis equal
xlabel('Fluor Y')
ylabel('Fluor C')

% Note that there is no time component in this data, so we cannot make
% cross correlation functions.
% Only the correlation coefficient can be measured.

%RmuY = corr(observableMu,observableY)
%RmuC = corr(observableMu,observableC)

% Slopes:
rMuY = corr(observableMu,observableY) % corr automatically normalizes by (std(observableMu)*std(observableY))
rMuC = corr(observableMu,observableC) % corr automatically normalizes by (std(observableMu)*std(observableC))

% And theoretically:
theoRMuY = noiseMuYCGainY / noiseMuYCGainMu
theoRMuC = noiseMuYCGainC / noiseMuYCGainMu

% Note that Noreen has sent me two websites with two nice explanations when
% the slope indeed is the same as the cov, and when not.
%
% - http://stats.stackexchange.com/questions/22718/what-is-the-difference-between-linear-regression-on-y-with-x-and-x-with-y/22721#22721
% - http://stats.stackexchange.com/questions/32464/how-does-the-correlation-coefficient-differ-from-regression-slope
