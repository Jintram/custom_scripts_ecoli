 % This script measures how well a probability estimate via KDE is by putting
 % errorbarbrs on the infered prob-distribution.
 % errorbars are created by bootstrapping from the measured data.
 
 % --------------------
 % Conclusions:
 % --------------------
 % works nicely.
 % The errorbars are rather small. 
 % For rare events they hit the zero-line (reasonable: from a n=1000 dataset
 %          we cannot conclude about the rare tails in prob distr p=10^-5 etc,
 %          so should be indistuingshable from =0)
 % robust to #bootstrapping redraws (low numbers as 100 &
 %          1000 gave very similar errorbars)
 % #gridN was little relevant
 
 
% Create "experimental" dataset
 datasize=1000; % #experimental datapoints
 
 % bimodal distribution:
 measuredvec=[randn(datasize/2,1);4+0.5*randn(datasize/2,1)]; % create "measured datapoints"
 
 
 %%
 
 numbootstrap=1000;    % # bootstrap repetitions
 % set grid parameters fix so that they don't vary for resampled datasets:
 rangedata=max(measuredvec)-min(measuredvec);
 minXgrid=min(measuredvec)-rangedata/10;    %same as default [for orig data,not necesarily resampled data]
 maxXgrid=max(measuredvec)+rangedata/10;
 gridN=2^12; % # datapoints in KDE grid 
 
 [~,probmeasured,gridX,~]=kde(measuredvec,gridN,minXgrid,maxXgrid); % get prob distr from experimental distr
  
 % a courageous step: create a large matrix with kde-prob for each
 % resampling
 % row1: prob1 of first resampling      prob1(x1)   prob1(x2) ...
 % row2: prob2 of 2nd resampling        prob2(x1)   prob2(x2) ...
 probboot_matrix=NaN(numbootstrap,gridN);
 
 % do resampling
 for run=1:numbootstrap
     % which idx are (re) drawn from orig distribution
     idxdraw=randi(datasize,datasize,1);
     bootvec=measuredvec(idxdraw);
     % kde on this redrawn vector
     [~,probboot,~,~]=kde(bootvec,gridN,minXgrid,maxXgrid);
     probboot_matrix(run,:)=probboot(:)';
 end
 
 
 
 % some plotting
 figure(1)
 clf 
 hold on
  for i=1:numbootstrap
     plot(gridX,probboot_matrix(i,:),'r')
  end
 plot(gridX,probmeasured,'k','LineWidth',3)
 plot(xlim,[0 0],'k')
 xlabel('x')
 ylabel('probdensity')
 
 
 
 % get errorbars (SEM) for some x-coordinates
 % SEM = stddev of all the calculated values (probs) with bootsrap)
 idxerr=[100:100:gridN];    %might need to adjust
 Xerr=gridX(idxerr);
 SEM=NaN(size(Xerr));
 for i=1:length(Xerr);
     semX=std(probboot_matrix(:,idxerr(i)))
     SEM(i)=semX;
 end
 
 
 
 figure(3)
 clf 
 hold on
errorbar(Xerr,probmeasured(idxerr),SEM)
plot(xlim,[0 0],'k')
  xlabel('x')
 ylabel('probdensity')