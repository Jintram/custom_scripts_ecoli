% Let's say there is 4 states: (A),(B),(C),(D).

% Calculate the error for the calculated probability for (A)
% (=standard error of the mean). use bootstrapping
% For (B),(C),(D) just repeat the steps below

% 1) Create a vector with length=sample size, (A) was measured -> 1, (sth
% else was measured) -> 0
% e.g. expdata=[1 1 0 0 0 0] means that have 6 data points, 2 times (A), 4
%                                    times sth else  ---->
%                                    mean(expdata)=0.33 is your probability
%                                    for which we want to get an error bar
%                                    so we need to resample the "mean"
expdata=zeros(1000,1); expdata(1:500)=1; % ADJUST!

% 2) Do bootstrapping (no need to adjust)
[bootresults,fullsampleidx]=bootstrp(10000,@mean,expdata);
% how to read: 10000 resamples. the "mean" of each resample will be in the
% output "bootresults", input was "expdata"
% fullsampleidx: not necessary for the moment

disp('first 3 resampled probabilities:')
bootresults(1:3)

 % 3) get the standard deviation of our resamplings = standard error of the
 % mean (no need to adjust)
 stderror=std(bootresults);
 disp(['your probability was: ' num2str(mean(expdata))]);
 disp(['standard error of the mean (errorbar) for the probability is: ' num2str(stderror)]);