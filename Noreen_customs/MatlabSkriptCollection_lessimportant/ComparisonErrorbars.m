%randn size: (N,1)
100               1000            10000             1000 uniform (rand

% theoretical mean
0               0               0                   0.5

%theoretical stddev
1               1               1                   0.3? or 0.28888?

%mean (sample mean)
-0.1250         -0.0516         0.0040              0.5030
% 4 replicates mean:
-0.1250   
-0.1108    
0.0689   
-0.1380
%their mean is: -0.0762

%std deviation
1.0647          0.9884          0.9889              0.2933
% stddev of the mean of the 4 replicates
0.0974

%SEM (std/sqrt(n))
0.1065          0.0313          0.0099              0.0093

%bootstrapping std dev of the means
%n=1000 draws
0.1063         0.0321           0.0101              0.0092
%n=10000 draws
0.1057         0.0312           0.0998              0.0093


%4  replicates of a uniform 100x1 distribution:
joint mean: 0.4991
std of the 4 means: 0.0220
SEM of 1st sample: 0.0271
bootstrap of 1st rep: 0.0273

%-> It is apparently equivalent (for the distribtuions investigated for which 
% each sample is representative (not systematically biased)
% to use Bootstrapping (std of the sampled means)
% or SEM (std/sqrt(n)) (n=sample size of 1 experiment)
% or Std of the means of several replicates
%but note that SEM underestimates (shown theoretically) for small n