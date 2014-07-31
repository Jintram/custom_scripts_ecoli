% Plots single cell traces of the production rate in dependence of
% cell cycle phase
% Additionally, the best fit of a (smoothed) step function (according 
% to the actually performed linear fit to fluorescence data) is added.
%
%
% PURPOSE OF THIS FUNCTION:
% choose e.g. three schnitzes which have a clear step in production rate
% and whose steps occurs at distinguishable time points
% These 3 traces will be color coded and the according fit will be plotted
% in the same color.
%
% NECESSARY ADJUSTMENTS
% 'schnitzStepDataVec': a vector that contains the schnitz number and the
%             step-fit for the production rate. Typically, this will be
%             'allstepsVec' from the code 'FitSmoothStepFunctionViaScanning.m'
%             It does not matter, if more schnitzes than used are written
%             into this vector.
%             content of one row: schnitz - phase (xjump) - prod (p0) - score -threshold
%             Note that p0 (from 'FitSmoot..') is normalized by colony average and will be
%             retransformed to absolute value within this function ('PlotSingle..').
% 'schnitznumrange': a vector that contains the schnitzes which will
%             actually be plotted
% 'schnitzUse': Schnitzcells structure which will be used. Typically
%             s_rm_fitTime (2012-05-08).
% 'NORMALIZE_P0': =1: all prodrate valus are normalized to the initial
%               fitted prodrate at phase=0. Reasonable, since absolute
%               values vary from cell-to-cell. Just dividing by average of the data points 
%               is not possible since the average is biased by the jump time).
%               =0: no normalization

% Further parameters (which typically will be the same for all plots) can
% be found below.
%
% ********** YOU NEED TO MAKE SURE THAT schnitzStepDataVec and schnitzUse REFER TO
% THE SAME COLONY!!! **********
%

% ------------------------------------------------------------
% ********* ADJUST **********
schnitzStepDataVec=allstepsVec;
schnitznumrange=[433 446 255];
% good schnitz candidates for 2012-05-08: early(433,182), middle(446,142), late(255)
NORMALIZE_P0=1;

% ** typically leave as is **
schnitzUse=s_rm_fitTime;
myfield='dY5';   % alternatively dC5, dY5_sum_dt etc
myphase='phase2_at_dY5_time';
ColorMatrix=[ 0 0 1; 0 0.8 0; 1 0 0];
% ------------------------------------------------------------

% add (random) colors if ColorMatrix is shorter than schnitznumrange
while size(ColorMatrix,1)<length(schnitznumrange)
    ColorMatrix=[ColorMatrix; rand rand rand];
end

% get total prodrate average to retransform initial production rates to
% tabsolute value (same calculation as in 'FitSmooth...')
yfpall=[];
for i=1:length(schnitzUse)
    if schnitzUse(i).useForPlot==1 & isempty(find(isnan(schnitzUse(i).dY5)))
        yfpall=[yfpall,schnitzUse(i).dY5];
    end
end
yfpmean=mean(yfpall);
% ********* RESCALE Prod0 ***********
schnitzStepDataVec(:,3)=schnitzStepDataVec(:,3)*yfpmean;
% *************************

% initialize figure
fig=figure(1);
clf
hold on
legendtext=[];

% -------------------------------
% loop over all schnitzes
% -------------------------------
for schnidx=1:length(schnitznumrange)
    schn=schnitznumrange(schnidx);
    
    % check if fit data over schnitz is present
    idx=find(schn==schnitzStepDataVec);
    if isempty(idx) | length(idx)>=2
        disp('Schnitz not found in schnitzDataStepVec or found more than once. Will stop here.');
    else
        % check if schnitz has complete cell cycle and useForPlot=1 (should
        % always be the case)
        if schnitzUse(schn).useForPlot==1 & schnitzUse(schn).completeCycle==1
            prodrate=schnitzUse(schn).(myfield);   %prodrate
            phase=schnitzUse(schn).(myphase); %phase
            
            % construct vector with fitted data
            prod0=schnitzStepDataVec(idx,3);
            jumpphase=schnitzStepDataVec(idx,2);            
            prodfit=zeros(size(phase));
            
            idxbefore=find(phase<jumpphase-0.001);  % before step  (always 0.01 buffer because of double variables (contrary to integer))
            if ~isempty(idxbefore)
                 prodfit(idxbefore)=prod0;
            end
            idxjump=find(phase<jumpphase+0.001 & phase>jumpphase-0.001);  % at step  (double variables! carefull with == ?!)
            if ~isempty(idxjump)
                 prodfit(idxjump)=prod0*1.5;
            end
            idxafter=find(phase>jumpphase+0.001);  % after step
            if ~isempty(idxafter)
                prodfit(idxafter)=prod0*2;
            end
            
            % normalization of prodrates
            if NORMALIZE_P0
                prodrate=prodrate/prod0;
                prodfit=prodfit/prod0;
            end
            
            % plotting
            figure(fig)
            g=plot(phase,prodrate,'.--','MarkerSize',15,'Color',ColorMatrix(schnidx,:));
            set(get(get(g,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            plot(phase,prodfit,'-','Color',ColorMatrix(schnidx,:),'LineWidth',2)
            
            % create legendtext
            if isempty(legendtext)
                legendtext={num2str(schn)};
            else
                legendtext=[legendtext, num2str(schn)];
            end
            
        end % end check complete CellCycle
     end % end check if schnitz exists in data vector
end % end loop over all schnitzes

figure(fig)
xlim([0 1])
legend(legendtext,'Location','NW')

clear idx idxbefore idxjump idxafter jumpphase prod0 phase prodrate prodfit1 
clear  myfield myphase ColorMatrix schn schnidx  g




