% This script is based on GetErrorbarsKDE_withBootstrap and is written specifically to obtain 
% errorboars for conditional distributions such as <YFP|mu> (in script: <field1|field3>
% (e.g. Fig.2 in paper)

% This script measures how well a probability estimate via KDE is by putting
 % errorbarbrs on the infered CONDITIONAL EXPECTATIONS
 % errorbars are created by bootstrapping from the measured data.
 
 % --------------------
 % Conclusions: BLUBB TO CHECK FOR NEW VERSION
 % --------------------
 % works nicely.
 % The errorbars are rather small. 
 % For rare events they hit the zero-line (reasonable: from a n=1000 dataset
 %          we cannot conclude about the rare tails in prob distr p=10^-5 etc,
 %          so should be indistuingshable from =0)
 % robust to #bootstrapping redraws (low numbers as 100 &
 %          1000 gave very similar errorbars)
 % #gridN was little relevant
 
 % ------------------------------
 % LOAD DATA
 % ------------------------------
 % load the experimental dataset from \\biofysicasrv\Users2\Walker\ExtrinsicNoise\Data_Collection\Data
 
 % **********  ADJUST ***********
 % choose fields of interest (column vectors!)
 % field names are according to NoiseDecomposition_viaKDE2.m
 field1=dC5_cycCor(:); % choose YFP or CFP rate (i.e. the field1or2, here always named field1)
 field3=mu_Third_cycCor(:); % always use growth rate
 % will calculate <field1|field3>
 
 cutoff_prob=0.001;   % min prob. for which the distribution is plotted (estimate gets bad for lower probs)
 numbootstrap=1000;    % # bootstrap repetitions
 % ******************************

 % normalize
 field1=field1/mean(field1);
 field3=field3/mean(field3);
 
% Create "experimental" dataset
datasize=length(field1); % #experimental datapoints
 
% % bimodal distribution:
% measuredvec=[randn(datasize/2,1);4+0.5*randn(datasize/2,1)]; % create "measured datapoints"
 
 
 %field1 -> X
 %field3 -> Y
 
 % set grid parameters fix so that they don't vary for resampled datasets (use default values):
 rangeXdata=max(field1)-min(field1);
 rangeYdata=max(field3)-min(field3);
 minXgrid=min(field1)-rangeXdata/4;
 maxXgrid=max(field1)+rangeXdata/4;
 minYgrid=min(field3)-rangeYdata/4;
 maxYgrid=max(field3)+rangeYdata/4;
 min_XY=[minXgrid minYgrid]
 max_XY=[maxXgrid maxYgrid]
 gridN=2^8; % # datapoints in KDE2d grid per axis
 
 %%
 % ---------------------------------
 % HARD CODE GRID PARAMETERS
 % ---------------------------------
 % This is a very messy step but since the option to fix the grid edges is
 % not implemented in NW_Probabilities_Grid_viaKDE2 it is the fastest
 % solution to hardcode the coordinates into kde2d_inclbandwidth.m between
 % line 103&104:
 % Enter: MAX_XY=max_XY; MIN_XY=min_XY (for the small letter variables
 % don't use the variable name but the actual number! DON'T FORGET TO
 % DELETE THESE LINES AFTER USAGE IN KDE2D_INCLBANDWIDTH!
 
 %%
 % ----------------------------
 % bootstrapping
 % ----------------------------
 
 % [bandwidth_f1f3,Probjoint_f1f3, meshX_f1f3, meshY_f1f3, Prob_f1,Prob_f3_fromf1, Prob_f1_given_f3, ...
 %   f1_grid, f3_grid, Mean_f1_given_f3, ~, ~, ~,increment_df3_fromf1]= ...
 %   NW_Probabilities_Grid_viaKDE2(field1,field3,gridN,1,1);
 
 % get conditional mean <field1|field3> from experimental distribution
 [~,~, ~, ~, ~,Prob_f3, ~, ~, f3_grid, Mean_f1_given_f3, ~, ~, ~,increment_df3]= ...
    NW_Probabilities_Grid_viaKDE2(field1,field3,gridN,1,1);
 % conditional prob distribution: Mean_f1_given_f3
 % growth rate axis: f3_grid
 

  
 % a courageous step: create a large matrix with conditional expectation for each
 % resampling
 % row1: cond.exp of first resampling      prob1(x1)   prob1(x2) ...
 % row2: cond.exp of 2nd resampling        prob2(x1)   prob2(x2) ...
 condexp_boot_matrix=NaN(numbootstrap,gridN);
 
 % do resampling
 for run=1:numbootstrap
     % which idx are (re) drawn from orig distribution
     idxdraw=randi(datasize,datasize,1);
     bootvec_field1=field1(idxdraw);
     bootvec_field3=field3(idxdraw);
     % kde on this redrawn vector
     
      [~,~, ~, ~, ~,~, ~, ~, ~, boot_Mean_f1_given_f3, ~, ~, ~,~]= ...
    NW_Probabilities_Grid_viaKDE2(bootvec_field1,bootvec_field3,gridN,1,1);
     
     condexp_boot_matrix(run,:)=boot_Mean_f1_given_f3(:)';
     
     if mod(run,20)==0
         disp(['Finished bootstrap repeat ' num2str(run)]);
     end
 end
 
 %%
 % -----------------------------------------
 % some plotting
 % ----------------------------------------
 %get idx together for which data will be plotted [minimum probability]
 idxplot=find(Prob_f3*increment_df3>cutoff_prob);
 
 figure(4)
 clf 
 hold on
  for i=1:numbootstrap
     plot(f3_grid(idxplot),condexp_boot_matrix(i,idxplot),'r')
  end
 plot(f3_grid(idxplot),Mean_f1_given_f3(idxplot),'k','LineWidth',3)
 plot(xlim,[0 0],'k')
 xlabel('x')
 ylabel('cond. exp.')
 
 
 
 % get errorbars (SEM) for some x-coordinates
 % !! SEM = stddev of all the calculated values (probs) with bootstrap !!
 idxerr=[6:10:2^8];    %might need to adjust
 % take overlap with idxplot
 idxerr=intersect(idxerr,idxplot)';
 
 Xerr=f3_grid(idxerr);
 SEM=NaN(size(Xerr));
 for i=1:length(Xerr);
     semX=std(condexp_boot_matrix(:,idxerr(i)));
     SEM(i)=semX;
 end
 
 
 
 figure(5)
 clf 
 hold on
 plot(f3_grid(idxplot),Mean_f1_given_f3(idxplot),'k','LineWidth',2)
errorbar(Xerr,Mean_f1_given_f3(idxerr),SEM,'k.','LineWidth',2)
%errorbar(Xerr,Mean_f1_given_f3(idxerr),SEM,'r.','LineWidth',1)
plot(xlim,[0 0],'k')
  xlabel('x')
 ylabel('cond. exp.')