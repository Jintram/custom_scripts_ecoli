% This skript is used in context with cell cycle paper (2013-02-06, NW)
% It compares the quality of different fits to cellular growth (length
% increase).
% For each cell with complete cell cycle (and useForPlot=1), an exponential
% fit, a linear fit and a bilinear fit (the latter via scanning) is
% attempted. The scores (deviationon from best fit) are stored and plotted
% EXACT DEFINITION 'SCORE' (updated 2013-11-13): 
% 1) Mean squared error, and the mean of the data (lengthvector) is normalized
%    to =1 (that is, the same criterion as in FitSmoothStepFunctionViaScanning
%    MSE = 1/n * Sum(deviations^2) / mean(data)^2  (Don't miss the last
%    square, otherwise normalized deviations individually without a square)
% 2) Akaike Information Criterion AIC. This criterion takes a penalty for #
%    of fitting paramters into account. For further info, check that
%    function.
%    Normalization of Data shouldn't matter (maybe if different schnitzes compared).
%    Here it's normalized to mean=1.
%    Lowest AIC -> Best model
%    Since only relative values of AIC matter, AIC(exponential fit) is set
%    to =0 and the relative differences are reported.
%
% General: LARGE SCORE = BAD FIT

% INPUT
% 'myschnitzcells'      which schnitzcells to use, typically s_rm_fitTime
% 'lambda_incr'          increment of linear growth rates for bilinear fit.
%                       Larger value -> fast calculation but inaccurate
% 'time_incr'           Increment between possible time points for switch 
%                       between the 2 linear rates in bilinear-fit    
% 'L0_incr'             Increment btw possible initial cell lengths for
%                       bilinear-fit
% 'PlotSchnitzes'       if not empty, plots fits and displays scores for
%                       these schnitzes
%
% OUTPUT
% 'fitScores_MSE'        contains schnitz numbers and MSE scores for different
%                      fits and further data on best bilinear fit
%                       schn1 --- score(exp) --- score(lin) --- score(bilin) --- lambda1 --- lambda2 --- L0 ---- timeswitch ---- phaseswitch
%                       schn2 --- score(exp) --- score(lin) --- score(bilin) --- lambda1 --- lambda2 --- L0 ---- timeswitch ---- phaseswitch
% 'fitScores_AIC'      same as fitScoresMSE, however with AIC scores
%
% Histogram with score distribution and average scores for each fit (ToDo??)
%
%
%
% ------------------------------------------
% NOTE!!!!!!!!!!!!! THE CURRENT DATA AND FIGURES ARE PROBABLY NOT UPDATED
% YET IN THE CELLCYCLE FOLDER!!!
% ------------------------------------------
%
% -----------------------------------------------------------------------------------

% ***********************************************************************************
% ******** ADJUST ****************
myschnitzcells=s_rm_fitTime;
%for bilinear fit of single traces:
lambda_incr=0.005/5;  %(0.02/5)    0.02: ca monolinear fit (linear growthr ate)
time_incr=2;            % (2)
L0_incr=2/50;    % (2/50)            2: ca monolinear fit (initial length)
PlotSchnitzes=[100:100:1000];%:2:160];
FORCELAMBDADOUBLE=0;  % =0 growth in first and second phase has independent rates. 
                      % =1 second rate is twice the first rate
DISPLAYBILINFITPROGRESSION=1; % display at which fitting step the bilin fit is (-> time estimate)
 
FITBINNEDDATA=1; % =1: only the binned data vector is fitted (at higher accuracy).
% For this, the vectors 'phasebins_forlength' and 'length_binned' have to be loaded
%for bilinear fit of binned data (.../CellCycle/Matlab/schnitzcells[...].m)
lambda_incr_binned=0.01; %(0.01)
time_incr_binned=0.01; %(0.01)
L0_incr_binned=0.005; % (should be 0.005)
% ***********************************************************************************
% ***********************************************************************************

% -----------------------------------------------------------------------------------
% ----------- initiate and set constants ---------------
% -----------------------------------------------------------------------------------
% initiate storage vector for quality of fits
%fitScores_MSE=zeros(0,4);
%fitScores_AIC=zeros(0,4);
fitScores_MSE=zeros(0,9);
fitScores_AIC=zeros(0,9);
% for AIC: regressionparameters:
K_lin=3;
K_exp=3;
K_bilin=5;

% ------------------------------------------------------------------------------------
% ----------- loop over all schnitzes (or a subset of interest) ------------
% If only the binend data is fitted, this loop is set to a trivial length
% of 1.
% ------------------------------------------------------------------------------------
schnitzloop=1:length(myschnitzcells); %adjust for debugging (blubb)
if FITBINNEDDATA
    schnitzloop=1;
end
for i=schnitzloop
    % let know where we are in loop
    if ~FITBINNEDDATA
        disp(['fitting schnitz ' num2str(i) '...'])
    else
        disp('fitting binned data ...')
    end
    clear birth division timevec lengthvec reltime continuoustime
    clear mu_exp L0_exp log_L0_exp lambda_lin L0_lin 
    clear fittedlength_exp fittedlength_lin fittedlength_bilin
    %load data (overwritten for binned fit)
    s=myschnitzcells(i);
    % fit data if (1) binned data or (2) completeCycle etc. of schnitz
    % exists
    if FITBINNEDDATA | (s.useForPlot==1 & s.completeCycle==1 & ~isnan(s.birthTime) & ~isnan(s.divTime)) %the latter 2 should always be true
        % ---------------------------------------------------------------------------------
        % --------- get data -----------
        % ---------------------------------------------------------------------------------
        if ~FITBINNEDDATA
            birth=s.birthTime;
            division=s.divTime;
            timevec=s.time;
            lengthvec=s.length_fitNew;
            reltime=timevec-birth;
            continuoustime=reltime(1):0.1:reltime(end);
        else
            % BINNEDDATA IS IN 'phase' NOT 'time' UNITS!
            birth=0; division=1;
            timevec=phasebins_forlength; %load from /CellCycle/Matlab/schnitzcells[...].m  %bins for the population-average length vector
            reltime=timevec;
            lengthvec=length_binned;  %load from /CellCycle/Matlab/schnitzcells[...].m  %population-averaged (binned) length vector
            continuoustime=0:0.01:1;
            lambda_incr=lambda_incr_binned;
            time_incr=time_incr_binned;
            L0_incr=L0_incr_binned;
            PlotSchnitzes=1; % set trivially to 1 (use same number as for trivial-loop-schnitz)
        end
        
                
        % ------------------------------------------------------------------------------------
        % ------ fit exponential -------
        % ------------------------------------------------------------------------------------
        loglengthvec=log2(lengthvec);  % linear fit to log2 of the length
        [fitparam_exp]=polyfit(reltime,loglengthvec,1);     % log2(L)=log2(L0)+mu*Time
        mu_exp=fitparam_exp(1); log_L0_exp=fitparam_exp(2);
        L0_exp=2^log_L0_exp;
        % fitted length values at each time point
        log_fitted_exp=polyval([mu_exp, log_L0_exp], reltime);
        fittedlength_exp=2.^(log_fitted_exp);
        % for plot
        continuousfittedlength_exp=2.^(polyval([mu_exp, log_L0_exp], continuoustime));
        %  ------------ deviation from fitted values -----------
        % MSE
        MSE_exp=sum((fittedlength_exp-lengthvec).^2); % sum of squared deviations
        MSE_exp=MSE_exp/(mean(lengthvec)^2); % normalize by cell length
        MSE_exp=MSE_exp/length(reltime); % normalize by number of data points (to get 'mean' instead of 'sum' SE)
        % AIC
        my_sigma2_est=MSE_exp; 
        my_n=length(lengthvec);
        my_K=K_exp;
        AIC_exp=AkaikeInformationCriterion(my_sigma2_est,my_n,my_K);
        
        % ------------------------------------------------------------------------------------
        % ------ fit linear -------
        % ------------------------------------------------------------------------------------
        [fitparam_lin]=polyfit(reltime,lengthvec,1);
        lambda_lin=fitparam_lin(1);
        L0_lin=fitparam_lin(2);
        % fitted length values at each time point
        fittedlength_lin=polyval([lambda_lin, L0_lin], reltime);
        % for plot
        continuousfittedlength_lin=polyval([lambda_lin, L0_lin], continuoustime);
        %  ------------ deviation from fitted values -----------
        % MSE
        MSE_lin=sum((fittedlength_lin-lengthvec).^2); % sum of squared deviations
        MSE_lin=MSE_lin/(mean(lengthvec)^2); % normalize by cell length
        MSE_lin=MSE_lin/length(reltime); % normalize by number of data points (to get 'mean' instead of 'sum' SE)
        % AIC
        my_sigma2_est=MSE_lin; 
        my_n=length(lengthvec);
        my_K=K_lin;
        AIC_lin=AkaikeInformationCriterion(my_sigma2_est,my_n,my_K);
        
        % ------------------------------------------------------------------------------------
        % ------ fit bilinear (via scanning)-------
        % ------------------------------------------------------------------------------------
        timestep=mean(diff(reltime)); % average timestep btw data points
        
        minL0_bilin=min(lengthvec)*0.9;  % heuristic: minimal length at birth time (before first data point!)
        maxL0_bilin=max(lengthvec);
        min_lambda1_bilin=min(diff(lengthvec))/(timestep);  % initial, lower lambda
        %min_lambda1_bilin=0;  % initial, lower lambda
        max_lambda1_bilin=max(diff(lengthvec))/(timestep); % lower than lambda2
        %  min_lambda2_bilin=0;  % sconde, higher lambda   % higher than lambda1
        max_lambda2_bilin=max(diff(lengthvec))/(timestep); % only relevant for flexible lambda2
        
        % remember increments: lambda_incr, time_incr, L0_incr
        
        bestdev=1000000; % use absurd high starting value for root mean squared deviation from fit
        best_piecewisepol=[];
        bestlambda1=[];    %fit parameter 1  
        bestlambda2=[];    %fit parameter 2  (if not set 2*lambda1)
        bestL0=[];     %fit parameter 3  
        besttimeswitch=[];     %fit parameter 4  
        bestphaseswitch=[];
        
       % % ****TEST *******
       % min_lambda1_bilin=0.5*lambda_lin;  % initial, lower lambda
       % max_lambda1_bilin=2*lambda_lin;
       % %  min_lambda2_bilin=0;  % sconde, higher lambda   % higher than lambda1
       % max_lambda2_bilin=2*lambda_lin;
       %  min_lambda1_bilin=min(diff(lengthvec))/(timestep);  % initial, lower lambda
        % ************
        
        
        % loop over all initial lengths L0
        for runL0=minL0_bilin:L0_incr:maxL0_bilin
            % create some output
            numL0=length(minL0_bilin:L0_incr:maxL0_bilin);
            numrunL0=round((runL0-minL0_bilin)/L0_incr);
            if DISPLAYBILINFITPROGRESSION
                disp(['bilinear fit. testing L0=' num2str(runL0) '. It''s ' num2str(numrunL0+1) ...
                  ' out of ' num2str(numL0) '.'])
            end
            % loop over all switch times for 2 linear subranges
            for runtimeswitch=reltime(1):time_incr:reltime(end) %(birth->0)
                % loop over all inital linear growth rates lambda1
                for lambda1run=min_lambda1_bilin:lambda_incr:max_lambda1_bilin
                    % flexible lambda2: loop over all second linear growth rates lambda2
                    % (larger than lambda1)
                    if ~FORCELAMBDADOUBLE
                        lambda2runrange=lambda1run:lambda_incr:max_lambda2_bilin;
                    else
                        % lambda2 is 2* lambda1 (fit constriction)
                        lambda2runrange=2*lambda1run; % one number
                    end
                    % loop over lambda2
                    for lambda2run=lambda2runrange
                        % create piecewise polynomial
                        L0switch=runL0+runtimeswitch*lambda1run; % length at switchtime
                        pptest=mkpp([0, runtimeswitch,reltime(end)], ...
                            [lambda1run,runL0; lambda2run, L0switch]);
                        % evaluate piecwise polynomial for measured time
                        % points and get quality of fit
                        currentfit=ppval(pptest,reltime);
                        currentdev=sum((currentfit-lengthvec).^2); % squared deviation
                        if currentdev<bestdev
                            bestdev=currentdev;
                            best_piecewisepol=pptest;
                            bestlambda1=lambda1run;
                            bestlambda2=lambda2run;
                            bestL0=runL0;
                            besttimeswitch=runtimeswitch;
                            bestphaseswitch=besttimeswitch/(division-birth);
                        end
                    end %loop over lambda2
                end %loop over lambda1
            end %loop over jump time
        end %loop over initial lanegth  ----> bilinear fit
        
        fittedlength_bilin=ppval(best_piecewisepol,reltime);
        continuousfittedlength_bilin=ppval(best_piecewisepol,continuoustime);
        %  ------------ deviation from fitted values -----------
        % MSE
        MSE_bilin=sum((fittedlength_bilin-lengthvec).^2); % sum of squared deviations
        MSE_bilin=MSE_bilin/(mean(lengthvec)^2); % normalize by cell length
        MSE_bilin=MSE_bilin/length(reltime); % normalize by number of data points (to get 'mean' instead of 'sum' SE)
        % AIC
        my_sigma2_est=MSE_bilin; 
        my_n=length(lengthvec);
        my_K=K_bilin;
        AIC_bilin=AkaikeInformationCriterion(my_sigma2_est,my_n,my_K);
        
        
        
        % -----------------------------------------------------------------------------
        % ---- only relative AICs matter ---------
        % -----------------------------------------------------------------------------
        AIC_lin=AIC_lin-AIC_exp;
        AIC_bilin=AIC_bilin-AIC_exp;
        AIC_exp=0;
        
        % -----------------------------------------------------------------------------
        % ----------- plot ---------
        % -----------------------------------------------------------------------------        
        
        % plot schnitz length including fits
        if ismember(i,PlotSchnitzes)
            figure
            clf
            hold on
            if FITBINNEDDATA
                title('Exp/Lin/Bilin Fitted binned length data. (MSE multiplied by 10^6)')
                xlabel('phase')
                ylabel('length [mum]')
            else
                title(['Exp/Lin/Bilin Fitted length data of schnitz ' num2str(i) '. (MSE multiplied by 10^6)'])
                xlabel('time after birth [min]') 
                ylabel('length [mum]')
            end
            plot(continuoustime,continuousfittedlength_exp,'-b','LineWidth',2)
            plot(continuoustime,continuousfittedlength_lin,'-r','LineWidth',2)
            plot(continuoustime,continuousfittedlength_bilin,'-g','LineWidth',2)
            plot(reltime,lengthvec,'.k','MarkerSize',15)
            %label1=['MSE_exp=' num2str(MSE_exp) ' AIC_exp=' num2str(AIC_exp)];
            %label2=['MSE_lin=' num2str(MSE_lin) ' AIC_lin=' num2str(AIC_lin)];    
            %label3=['MSE_bilin=' num2str(MSE_bilin) ' AIC_bilin=' num2str(AIC_bilin)];    
            label1=sprintf('MSE_exp  = %4.3f      AIC_exp = %4.3f',MSE_exp*10^6,AIC_exp);
            label2=sprintf('MSE_lin   = %4.3f    AIC_lin = %4.3f',MSE_lin*10^6,AIC_lin);
            label3=sprintf('MSE_bilin = %4.3f    AIC_bilin = %4.3f',MSE_bilin*10^6,AIC_bilin);
            le=legend(label1,label2,label3, 'Location','NW');
            set(le,'Interpreter','None');
        end
            
        
        fitScores_MSE=[fitScores_MSE; i, MSE_exp, MSE_lin, MSE_bilin ...
            bestlambda1, bestlambda2, bestL0, besttimeswitch, bestphaseswitch];
        fitScores_AIC=[fitScores_AIC; i, AIC_exp, AIC_lin, AIC_bilin ...
            bestlambda1, bestlambda2, bestL0, besttimeswitch, bestphaseswitch];
    end % if schnitz has complete cell cycle ...
    
end % loop over all schnitzes