function stepschnitzes = FindSteps_cellCycle(schnitzarray,fieldphase,fieldcolor,varargin)
% find schnitzes which seem to have a step in the cell cycle dependence of
% fieldcolor (typically yfp production rate
% schnitzarray, fieldphase & fieldcolor are obtained by
% GetCellCycleData_for_PhasePlot.m
% schnitznumber range can be given optionally (via varargin), e.g. [400
% 420]
%
% Currently used conditions are:
% (1) significant increase: smooth data. mean of first 2 datapts must be at
%     least XX* smaller than mean of last 2 datapts (XX is adjustable)
% (2) no decrease: smooth data. max YY datapts may exist (noise) for which subsequent one
%     has a smaller value (YY is adjustable)
% (3) (sudden) step: it must happen within 3 datapoints (2 time steps).
%     difference between one and third datapoint must be at least
%     ZZ*mean(fieldcolor) (and of course positive) (ZZ is adjustable)


% get manual schnitzes range if given
manrange=0;
if length(varargin)>0
    schnitzrange=varargin{1};
    schnitzmin=schnitzrange(1);
    schnitzmax=schnitzrange(2);
    manrange=1;
end
    

% ***** FINETUNE **********
Threshold_SignIncrease=1.5;     % how much must fieldcolor increase
Threshold_DecreasePoints=1;     % how many datapoints may decrease in value
Threshold_StepSize=0.8;         % how large is step (compared to mean(fieldcolor)
% *************************

% ***** PLOTOPTIONS & Further Options (Debugging) *******
NumPlotsMax=20; % how many schnitzes are plotted? if =0, no plotting
PlotAllSchnitzes=1; %plots all schnitzes, no matter whether step determined or not
                    % if =1, NumPlotsMax refers to total # of plots (and
                    % not step-plots)
DISPLAYDATA=1;      % print data of subfunctions (debugging)
% *************************


% find unique schnitznumbers
sch=unique(schnitzarray);
% if manual schnitzes range, get subselection
if manrange
    idxman=find(sch>=schnitzmin & sch<=schnitzmax);
    sch=sch(idxman);
end
% array for schnitzes with step
stepschnitzes=[];

% loop over all schnitzes
for i=1:length(sch);
    schnitzrun=sch(i);
    idx=find(schnitzarray==schnitzrun);
    currphase=fieldphase(idx);
    currcolor=fieldcolor(idx);
    % some datapoints exist (min length 3)
    if length(currphase)>=3
        % (1) increase significant
        yesno1=SubFct_SignificantIncrease(currcolor,Threshold_SignIncrease,DISPLAYDATA);
        % (2) few datapoints (1) with decrease
        yesno2=SubFct_LimitDecrease(currcolor,Threshold_DecreasePoints,DISPLAYDATA);
        % (3) large stepsize
        yesno3=SubFct_LargeStepSize(currcolor,Threshold_StepSize,DISPLAYDATA);
        
        % if all conditions are met, add schnitz number to array of
        % stepschnitzes
        if yesno1 & yesno2 & yesno3
            stepschnitzes=[stepschnitzes, schnitzrun];
        end
        
        if DISPLAYDATA
            disp(['schn=' num2str(schnitzrun) '    yesno1=' num2str(yesno1) ...
                ' yesno2=' num2str(yesno2) ' yesno3=' num2str(yesno3)])
            disp(' ')
        end
        
    end
        
    
end

% --------- PLOTTING -----------
% only plot stepschnitzes
if ~PlotAllSchnitzes
    numplots=min(NumPlotsMax,length(stepschnitzes)); % number of figures
    if NumPlotsMax>0 & length(stepschnitzes)>NumPlotsMax
        disp(['more step-schnitzes found than plots chosen. Will display first ' num2str(NumPlotsMax) ' found schnitzes']);
    end
    for k=1:numplots
        schnpl=stepschnitzes(k);
        figure
        idx=find(schnitzarray==schnpl);
        plot(fieldphase(idx),(fieldcolor(idx)),'.-r')
        %hold on; plot(ph(idx),(gfp(idx)),'.-g')
        title([num2str(schnpl)])
    end
else
% plot all schnitzes
    numplots=min(NumPlotsMax,length(sch)); % number of figures
    if NumPlotsMax>0 & length(sch)>NumPlotsMax
        disp(['more schnitzes examined than plots chosen. Will display first ' num2str(NumPlotsMax) ' schnitzes.']);
    end
    for k=1:numplots
        schnpl=sch(k);
        figure
        idx=find(schnitzarray==schnpl);
        idxstep=find(stepschnitzes==schnpl);
        if ~isempty(idxstep)
            plot(fieldphase(idx),(fieldcolor(idx)),'.-r')
            title([num2str(schnpl) ' step'])
        else
            plot(fieldphase(idx),(fieldcolor(idx)),'.-b')
            title([num2str(schnpl) ' no'])
        end
    end
end % end plotting

end % end main function

% ********* TEST FUNCTIONS FOR STEP *************    
% yesno=1 -> test passed
% tesno=0 -> test failed
% ***********************************************

    % -------------------------------------
    % ---- (1) checks if significant increase btw first and last datapts
    function yesno=SubFct_SignificantIncrease(field,MinIncrease,DISPLAYDATA)
        field=smooth(field);
        initvalue=mean(field(1:2));
        endvalue=mean(field(end-1:end));
        frac=endvalue/initvalue;
        % yes=increase significant  ->1
        % no=low increase -> 0
        yesno=(frac>=MinIncrease);
    end

    % -------------------------------------
    % ---- (2) checks that max XX (e.g. 1) datapts decrease
    function yesno=SubFct_LimitDecrease(field,MaxNumber,DISPLAYDATA)
        field=smooth(field); %blubb
        fielddiff=diff(field);
        % *** blubb
        totaldecr=sum(fielddiff(fielddiff<0));
        reldecr=totaldecr/mean(field);
        totalincr=sum(fielddiff(fielddiff>0));
        numdecr=length(find(fielddiff<0));
        if DISPLAYDATA
            disp(['totaldecrease='  num2str(totaldecr) ...
                '  rel. decr =' num2str(reldecr) '   # decr =' ...
                num2str(numdecr) '   decr/incr=' num2str(totaldecr/totalincr)])
        end
        % ***
        decreasepts=find(fielddiff<0);
        %if isempty(decreasepts) | length(decreasepts)<=MaxNumber
        if abs(reldecr)<0.4
            yesno=1;
        else
            yesno=0;
        end
    end

    % -------------------------------------
    % ---- (3) checks if step size is large enough
    function yesno=SubFct_LargeStepSize(field,MinStepSize,DISPLAYDATA)
        meanfield=mean(field);
        yesno=0;
      
        for j=1:length(field)-2;
            before=field(j);
            after=field(j+2);
            incr=after-before;
            relincr=incr/meanfield;
            if relincr>MinStepSize;
                yesno=1;
                break
            end
        end
    end


        
        
        