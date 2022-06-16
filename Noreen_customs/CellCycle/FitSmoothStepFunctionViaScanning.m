% dumb method to find best continuous step function to single cell
% produciton traces
%-> scan parameter range and search for least deviation 
% 2 parameters:
% p0=initial production rate (doubles at the end)
%   - restriciton used: prod rate has to be 2*p0 at the end of cell cycle
% xjump=timepoint at which production rate doubles (=x in script)
%   - simplifications used: jump can only occur at a datapoint (ca 7
%      time/phase points) and not in between.
%      This implies that p (prod rate) at this time point is 1.5*p0 and
%      2*p0 at the subsequent one (if intermediate jumps are allowed, 2
%      instead of one datapoint may show intermediate p values
%
% NOTES: pages 164ff
%
% scan range:
% p0: scan from minimal to maximal p of the specific trace. stepsize is
%   given manually (maybe for further improvement make more detailed scanning
%   in 'good' region
% xjump: first till last phase time point 
% ******************************* 
% I think it is not necessary to use parent and daughter info
% as long as the jump is restricted to actual datatimepoints?   ???? 
% ********************************
%
% ***********************
% AT THE MOMENT THE (SQUARE) DEVIATION/SCORE FROM THE FIT IS NORMALIZED BY
% THE SQUARE OF THE AVERAGE PROD RATE OF THAT INDIVIDUAL TRACE
% FURTHER NORMALIZATION BY THE NUMBER OF DATAPOINTS (not squared)
% ---- SHOULD IT BE SQUARED???? -------
% ***********************
%
% Shape of fitting function
%
%    2p0|             *  *  *  *
%       |          *
%     p0| *  *  *  |
%       |          |
%       |__________|____________________
%                  xjump
%
%
% *** OUTPUT ***
%  'allstepsVec'  containing: schnitz - phase (xjump) - prod (p0) - score -threshold - timejump [min] - interdivtime [min]
%  -threshold - time[min] of jump (rel. to birth) - lifetime of cell [min]
%  'stepFraction' : fraction of schnitzcells which is considered to have a
%                   step function trace (dependent on threshold)
%
% -----------------------------------------------------------------------

% **** ADJUST ****
pincr=0.02; %increment of production rates tested (usually prod rate is normalized to average=1)
PLOTFIT=0;
PLOT2ND3RDFIT=0;
PLOTONLYSTEPSCHNITZES=1;
SAVESTEPTRACEFIGURES=XXX1; % the other figures can currently not be saved automatically
SAVEDIR='XXX\\storage01\data\AMOLF\users\walker\CellCycle\Data_and_Analysis\StepTrace_Schnitzes_and_Distribution\SomData_UseFullTimeRange\AllIndividualStepTraces_incl_Fit\';
schnitzUse=s_rm_fitTime;
schnitznumrange=1:length(schnitzUse); % earlysteps';
ExcludeFirstStep=0;  % if=1, the step xjump must not occur at first datapoint (makes sense since typical dbl time =100 min, 
                     %  >=6 frames per cell -> first data point at <20min.
                     %  cell cannot yet have doubled genome at yfp position
                     %  + maturation time (only exception simultaneous
                     %  multiple replication forks at fast growth)
MinMinutes=0;        % minimal amount of minutes after division which is 
                     % admitted as gene dbl time (somewhat redundant to ExcludeFirstStep).
                     % Theoretically, the first ca 25min there is no
                     % YFP-dbl possible (single replication fork)
                     % =0: allow all step times
DeviationThreshold=0.02; % until which deviation is the data still considered a step function?
% ****************


% get total average and later normalize by it;
yfpall=[];
for i=1:length(schnitzUse)
    if schnitzUse(i).useForPlot==1 & isempty(find(isnan(schnitzUse(i).dY5)))
        yfpall=[yfpall,schnitzUse(i).dY5];
    end
end
yfpmean=mean(yfpall);
% prepare vector that stores all stepschnitzes including steptime, p0 and score
allstepsVec=zeros(0,7);   % schnitz - phase (xjump) - prod (p0) - score -threshold

% initialize to calculate fraction of step traces
stepSchnitzes=0;
investigatedSchnitzes=0;

% ----------------------------------------------------------
%  LOOP OVER ALL SCHNITZES 
% ----------------------------------------------------------

for schnitznum=schnitznumrange
    %schnitznum=440; %debugging
    ss=schnitzUse(schnitznum);
    % complete cycle?
    if ss.completeCycle==1 & ss.useForPlot==1

        %ph_schnitz=[1 2 3 4 5 6 7]; %debugging
        %prod_schnitz=[1 0.8 1.2 1 1.5 2 2]; %debugging
        ph_schnitz=ss.phase2_at_dY5_time;
        prod_schnitz=ss.dY5/yfpmean;
        time_schnitz=ss.dY5_time;
        birth_time_schnitz=ss.birthTime;
        div_time_schnitz=ss.divTime;
        life_time_schnitz=div_time_schnitz-birth_time_schnitz;

        % get range of steps that production rate can have for fit function
        pmin=min(prod_schnitz);
        pmax=max(prod_schnitz);
        psteps=0.7*pmin:pincr:0.5*1.3*pmax;  % initial rate is at most half of the max rate ever measured
                                % little bit arbitrary range

        numph=length(ph_schnitz);
        numprod=length(psteps);

        % vector that saves all fit parameters and deviations/scores
        if ~ExcludeFirstStep
            deviation=zeros(numph*numprod,3);  % phase (xjump) - prod p0 - deviation    ........per row
        else
            deviation=zeros((numph-1)*numprod,3); 
        end

        % ----------------------------------------------------
        % loop over all jump time points 
        % ----------------------------------------------------
        if ~ExcludeFirstStep
            phaserangeidx=1:length(ph_schnitz);
        else
            phaserangeidx=2:length(ph_schnitz);
        end
        for phrunidx=phaserangeidx
            phrun=ph_schnitz(phrunidx);
            % ----------------------------------------------------------
            % loop over all initial production rates
            % ----------------------------------------------------------
            for prodrunidx=1:numprod
                prodrun=psteps(prodrunidx);

                % calculated values with current fit function
                prod_calc=zeros(1,length(ph_schnitz));
                idx=find(ph_schnitz<phrun);  % before step
                if ~isempty(idx)
                    prod_calc(idx)=prodrun;
                end
                idx=find(ph_schnitz<phrun+0.01 & ph_schnitz>phrun-0.01);  % at step  (double variables! carefull with == ?!)
                if ~isempty(idx)
                    prod_calc(idx)=prodrun*1.5;
                end
                idx=find(ph_schnitz>phrun);  % after step
                if ~isempty(idx)
                    prod_calc(idx)=prodrun*2;
                end

                %calculate deviation from experimental data
                currentdev=(prod_calc-prod_schnitz);
                currentdev=sum(currentdev.^2);
                % ******** NORMALIZE BY TRACE AVERAGE PROD RATE^2 and # OF DATAPOINTS **********
                currentdev=currentdev/(mean(prod_schnitz)^2*length(prod_schnitz));
                %currentdev=currentdev/(mean(prod_schnitz)^2*(length(prod_schnitz))^2);

                % where to write into deviation score matrix
                if ~ExcludeFirstStep
                    deviationidx=(phrunidx-1)*numprod+prodrunidx;
                else
                    deviationidx=(phrunidx-2)*numprod+prodrunidx;
                end
                deviation(deviationidx,:)=[ phrun, prodrun, currentdev];
            end
        end  % end loop over all parameters

        % --------------------------------------------
        % EXTRACT BEST FIT
        % --------------------------------------------
        deviationsorted=sortrows(deviation,3);
        fit1=deviationsorted(1,:);  % best fit  % phase - p0 - score
        
        % if minimal delay after division is set, search for first fit
        % satisfying this condition
        if MinMinutes>0  %cutoff
            birthTime=ss.birthTime;
            reltime_schnitz=ss.dY5_time-birthTime;
            idxtime=find(reltime_schnitz>=MinMinutes);
            minPhase=min(ph_schnitz(idxtime));
            fitidx=1; % just necessary in case fit is allright -> for calc of fit2
            if fit1(1)<minPhase
                for fitidx=2:length(deviationsorted)
                    fit1=deviationsorted(fitidx,:);
                    if fit1(1)>=minPhase;
                        break
                    end
                end
            end    
        end
        if fit1(3)<=DeviationThreshold
            ISSTEP=1;
            % get time of jump
            phasejump=fit1(1);
            idxtime2=find(ph_schnitz>fit1(1)-0.01 & ph_schnitz<fit1(1)+0.01); % dbl 
            timejump=time_schnitz(idxtime2)-birth_time_schnitz;
            clear phasejump idxtime2
            allstepsVec=[allstepsVec; schnitznum, fit1(1), fit1(2), fit1(3), DeviationThreshold, timejump, life_time_schnitz];
        else
            ISSTEP=0;
        end
        
        if PLOTONLYSTEPSCHNITZES
            PLOTSTEPAPPROVED=ISSTEP;
        else PLOTSTEPAPPROVED=1;
        end
        
        investigatedSchnitzes=investigatedSchnitzes+1;
        if ISSTEP
            stepSchnitzes=stepSchnitzes+1;
        end

        % ------------------------------
        % plotting
        % ------------------------------
        if PLOTFIT & PLOTSTEPAPPROVED
            prodfit1=zeros(size(ph_schnitz));
            idx=find(ph_schnitz<fit1(1));  % before step
            if ~isempty(idx)
                 prodfit1(idx)=fit1(2);
            end
                 idx=find(ph_schnitz<fit1(1)+0.01 & ph_schnitz>fit1(1)-0.01);  % at step  (double variables! carefull with == ?!)
            if ~isempty(idx)
                 prodfit1(idx)=fit1(2)*1.5;
            end
                idx=find(ph_schnitz>fit1(1));  % after step
            if ~isempty(idx)
                prodfit1(idx)=fit1(2)*2;
            end
        end

        if PLOT2ND3RDFIT & PLOTSTEPAPPROVED
            fit2=deviationsorted(2,:);  % 2nd fit
            if MinMinutes>0  % cutoff
                fit2=deviationsorted(fitidx+1,:);  % start one after fit1
                fitidx2=fitidx+1; % just necessary in case fit2 is allright -> for calc of fit3
                if fit2(1)<minPhase
                    for fitidx2=fitidx+1:length(deviationsorted)
                        fit2=deviationsorted(fitidx2,:);
                        if fit2(1)>=minPhase;
                            break
                        end
                    end
                end 
            end
            
            prodfit2=zeros(size(ph_schnitz));
            idx=find(ph_schnitz<fit2(1));  % before step
            if ~isempty(idx)
              prodfit2(idx)=fit2(2);
            end
            idx=find(ph_schnitz<fit2(1)+0.01 & ph_schnitz>fit2(1)-0.01);  % at step  (double variables! carefull with == ?!)
            if ~isempty(idx)
                prodfit2(idx)=fit2(2)*1.5;
            end
            idx=find(ph_schnitz>fit2(1));  % after step
            if ~isempty(idx)
                prodfit2(idx)=fit2(2)*2;
            end

            fit3=deviationsorted(3,:);  % 3rd fit
            if MinMinutes>0  % cutoff
                fit3=deviationsorted(fitidx2+1,:);  % start one after fit1
                if fit3(1)<minPhase
                    for fitidx3=fitidx2+1:length(deviationsorted)
                        fit3=deviationsorted(fitidx3,:);
                        if fit3(1)>=minPhase;
                            break
                        end
                    end
                end 
            end
            prodfit3=zeros(size(ph_schnitz));
            idx=find(ph_schnitz<fit3(1));  % before step
            if ~isempty(idx)
                prodfit3(idx)=fit3(2);
            end
            idx=find(ph_schnitz<fit3(1)+0.01 & ph_schnitz>fit3(1)-0.01);  % at step  (double variables! carefull with == ?!)
            if ~isempty(idx)
                prodfit3(idx)=fit3(2)*1.5;
            end
            idx=find(ph_schnitz>fit3(1));  % after step
            if ~isempty(idx)
                prodfit3(idx)=fit3(2)*2;
            end
        end  % end PLOT2ND3RDFIT

        if PLOTFIT & PLOTSTEPAPPROVED
            figure
            clf
            hold on
            title(['schnitz ' num2str(schnitznum) '.   squared&normed scores in legend']);
            if ~PLOT2ND3RDFIT
                plot(ph_schnitz,prodfit1,'-r','LineWidth',2)
                legend([num2str(fit1(3))])
            else
                %plot(ph_schnitz,prodfit1,'-r','LineWidth',2)
                plot(ph_schnitz,prodfit1,'.-r','LineWidth',2,'MarkerSize',15)
                plot(ph_schnitz,prodfit2,'-','Color',[1 0.5 0])
                plot(ph_schnitz,prodfit3,'-','Color',[1 0.8 0])
                legend([num2str(fit1(3))], [num2str(fit3(3))],[num2str(fit3(3))],'location','NW') 
            end
            plot(ph_schnitz,prod_schnitz,'.-b','MarkerSize',20,'LineWidth',2)
            xlim([0 1])
            
            %BLUBB start 2016-01-04 Quick adjustments to make the reviewers
            %happy. Check single cell step traces
            title('')
            set(gcf,'OuterPosition',[100 100 300 300])
            xlabel('phase')
            ylabel('prod. rate')
            %if ISSTEP
            %    annotation(gcf,'textbox',[0.77 0.14 0.12 0.06],'String',{'Step'});
            %end
            legend off
            if SAVESTEPTRACEFIGURES
                % save the figures
                figname=['SingleTrace_Schnitz' str3(schnitznum)];
                savefig(gcf,[SAVEDIR figname])
                saveSameSize(gcf,'format','png','file',[SAVEDIR figname '.png']) %.pdf does not work with same size
                saveSameSize(gcf,'format','eps','file',[SAVEDIR figname '.eps'])
            end
            %BLUBB end
            
            % BLUBB
           % if ~ISSTEP & fit1(3)<=0.03
           %     annotation(gcf,'textbox',[0.77 0.14 0.12 0.06],'String',{'Step0.03'});
           % end

        end % end PLOTFIT

        % ----------------------- end plotting --------------
        
    end % end completeCycle
end % end loop over all schnitzes

% -------------------------------------------------
% PLOT HISTOGRAM AND CUMULATIVE FCT OF STEPS
% -------------------------------------------------

figure
clf
%hist(allstepsVec(:,2),[0.03333:0.0666:1]) % first bar should not be centered around 0, but have it as left corner
hist(allstepsVec(:,2),[0.05:0.1:1])
hold on
xlim([0 1])
title('jump times')

% output: fraciton of stepschnitzes
disp('-----------------------------')
stepFraction=stepSchnitzes/investigatedSchnitzes;
disp(['Investigated Schnitzes: ' num2str(investigatedSchnitzes) ... 
    '.  Max tolerated deviation: ' num2str(DeviationThreshold)]);
disp(['Fraction considered as step-traces: ' num2str(stepFraction) ...
     ' (abs number: ' num2str(stepSchnitzes) ')']);
disp('-----------------------------')

% ********************************
% CALCULATE RECONSTITUTED CELL CYCLE DEPENDENCE OF PRODUCTION RATE
% *********************************
% - uses jumptime distribution of allstepsVec(:,2)
% - works with actual steps/jumps and not linear increase
%
% for initial production rate p0, three possibilities are used:
% 1) all start off with p0=1 (different normalization!)
% 2) p0 is chosen such that for each jumptime the total produced
%     Fluorescence protein =1 (use actual jump, not fitted lin increase)
% 3) use experimental p0 distribution, i.e. allstepsVec(:,3)


jumpxVec=allstepsVec(:,2);     % vector with all jump times
p0_case1=ones(size(jumpxVec)); % 1) 
p0_case2=1./(2-jumpxVec);      % 2) (see notes p166)
p0_case3=allstepsVec(:,3);

overallProd_case1=zeros(200,1);
overallProd_case2=zeros(200,1);
overallProd_case3=zeros(200,1);
overallPhase=1/200:1/200:1;

for i=1:length(overallProd_case1)
    idxbeforejump=find(jumpxVec>overallPhase(i));
    idxafterjump=find(jumpxVec<=overallPhase(i));
    overallProd_case1(i)=sum(p0_case1(idxbeforejump))+2*sum(p0_case1(idxafterjump));
    overallProd_case2(i)=sum(p0_case2(idxbeforejump))+2*sum(p0_case2(idxafterjump));
    overallProd_case3(i)=sum(p0_case3(idxbeforejump))+2*sum(p0_case3(idxafterjump));
end
overallProd_case1=overallProd_case1/mean(overallProd_case1);
overallProd_case2=overallProd_case2/mean(overallProd_case2);
overallProd_case3=overallProd_case3/mean(overallProd_case3);
    
figure
clf
hold on
plot(overallPhase,overallProd_case1,'b','LineWidth',2)
plot(overallPhase,overallProd_case2,'r','LineWidth',2)
plot(overallPhase,overallProd_case3,'Color',[0 0.8 0],'LineWidth',2)  %green
legend('same initial','integral 1','experim','Location','NW')
title(['StepFraction: ' num2str(stepFraction) ' (abs number: ' num2str(stepSchnitzes) ')'])
