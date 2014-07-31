% ----------------------------------------------------------------------
% ------- BOTH SIGNALS AT SAME TIME (DIFFERENT TIME: SEE NEXT CELL)
% ----------------------------------------------------------------------

mup_vec1=[];  % (mu, prodrate)
mup_vec2=[];
Delta_mup_vec=[];
muTimeField='Y_time';
rateTimeField='dY5_time';
rateField='dY5_cycCor';
% -----------------------------------
% Loop over all schnitzcells
% -----------------------------------
for i=1:length(testschnitzcells)
    s=testschnitzcells(i);
    % ----------------------------------------
    % Loop over all time points of specific schnitz
    % ---------------------------------------
    for idx1=1:length(s.(rateTimeField))
        mutimeidx1=find(s.(muTimeField)==s.(rateTimeField)(idx1));
        p1=s.(rateField)(idx1);
        mu1=s.muP15_fitNew_cycCor(mutimeidx1);
        % ------------------------------------
        % if not last entry, schnitz has also data for next time point 
        % -------------------------------------
        if idx1<=length(s.(rateTimeField))-1
            idx2=idx1+1;
            mutimeidx2=find(s.(muTimeField)==s.(rateTimeField)(idx2));
            if ~isempty(mutimeidx1) & ~isempty(mutimeidx2)  % should always be true (mu contains more data than the rates)
                p2=s.(rateField)(idx2);
                mu2=s.muP15_fitNew_cycCor(mutimeidx2);
                if ~isnan(p1+p2+mu1+mu2) % 'NaN' is true if any summand is NaN
                    mup_vec1=[mup_vec1; mu1, p1];
                    mup_vec2=[mup_vec2; mu2, p2];
                    Delta_mup_vec=[Delta_mup_vec; mu2-mu1, p2-p1];
                end
            end
        % ---------------------------------------
        % for last time point, take data from daugthers
        % ---------------------------------------
        else
            % check if daughters exist
            if s.D~=0 & s.E~=0
                d1=testschnitzcells(s.D);
                d2=testschnitzcells(s.E);
                % check if daughters have data (time vector min length = 1)
                if length(d1.(rateTimeField))>0 & length(d2.(rateTimeField))>0 & d1.(rateTimeField)(1)==d2.(rateTimeField)(1)
                    mutimeidx2d1=find(d1.(muTimeField)==d1.(rateTimeField)(1));
                    mutimeidx2d2=find(d2.(muTimeField)==d1.(rateTimeField)(1)); %both should be the same
                    if ~isempty(mutimeidx2d1) & ~isempty(mutimeidx2d2) & mutimeidx2d1==mutimeidx2d2
                        p2=d1.(rateField)(1)+d2.(rateField)(1); % sum up rates of daughters
                        mu2=d1.muP15_fitNew_cycCor(mutimeidx2d1)*d1.birthLengthRatio15+ ... 
                            d2.muP15_fitNew_cycCor(mutimeidx2d2)*d2.birthLengthRatio15; %average growth rates
                        if ~isnan(p1+p2+mu1+mu2) % 'NaN' is true if any summand is NaN
                            mup_vec1=[mup_vec1; mu1, p1];
                            mup_vec2=[mup_vec2; mu2, p2];
                            Delta_mup_vec=[Delta_mup_vec; mu2-mu1, p2-p1];
                        end
                    end
                end
            end
        end
        % --- end of: if i<=length(s.(rateTimeField))-1
    end
    % --- end of: loop over time points
end
% end of: loop over schnitzes
                    
            
%%            
                
% ----------------------------------------------------------------------
% ------- SIGNALS AT DIFFERENT TIME (DOES NOT AT DELTA_MU, DELTA_P)
% ---- SEARCHES (so far) ONLY ONE DAUGHTER/PARENT GENERATION, NOT MORE!!!
% ----------------------------------------------------------------------

% ********** ADJUST *************************************
muTimeField='G_time';
rateTimeField='dG5_time';
rateField='dG5_cycCor';
fitTime=[200 600];  % condition on growth AND production rate
timeShiftNum=0; % number of time steps that prodrate-signals shall differ from growth rate signal
    % =0: signals taken at same time
    % <0: growth rate taken at earlier time (=-1 one fluotimeframe earlier)
    % >0: growth rate taken at later time
    % e.g. for ribosomeL31 case: take ca. =4

testschnitzcells=s_rm; %rm schnitzes...
% ********************************************************

mup_vec1=[];  % (mu, prodrate)

%get time difference in minutes (allow tolerance of XX minutes)
%takes specific schnitz!
for i=1:length(testschnitzcells)
    s=testschnitzcells(i);
    if length(s.(rateTimeField))>=2
        DeltaTime=s.(rateTimeField)(2)-s.(rateTimeField)(1);
        break
    end
end
timeShift=timeShiftNum*DeltaTime; %time shift of signals in [min]
timeTolerance=0.5;
if timeShiftNum>=0 % for later distinction whether parent or daughters to look for
    BoolMuLate=1;
else
    BoolMuLate=0;
end

% -----------------------------------
% Loop over all schnitzcells
% -----------------------------------
for i=1:length(testschnitzcells) %blubb
    s=testschnitzcells(i);
    if (~existfield(s,'useForPlot') | s.useForPlot==1)
        % ----------------------------------------
        % Loop over all time points of specific schnitz
        % ---------------------------------------
        for idx1=1:length(s.(rateTimeField))
            ProdTime=s.(rateTimeField)(idx1); % time point of production rate
            MuTime=ProdTime+timeShift;
            % ----------------------------------------
            % Are mu time and rate time within fitTime?
            % ---------------------------------------
            if fitTime(1)<=ProdTime & fitTime(1)<=MuTime & fitTime(2)>=ProdTime & fitTime(2)>=MuTime               
                p1=s.(rateField)(idx1);
                % ------------------------------------
                % If Schnitz has data at 'MuTime'
                % -------------------------------------
                mutimeidx1=find(s.(muTimeField)>MuTime-timeTolerance & s.(muTimeField)<MuTime+timeTolerance);
                if ~isempty(mutimeidx1)
                    mu1=s.muP15_fitNew_cycCor(mutimeidx1);
                    if ~isnan(p1+mu1)
                        mup_vec1=[mup_vec1; mu1, p1, i];
                    end
                % ---------------------------------------
                % Elseif: Schnitz has no data of Mu. Mu is later, search daughters
                % ---------------------------------------    
                else
                    if BoolMuLate
                        % check if daughters exist
                        if s.D~=0 & s.E~=0
                            d1=testschnitzcells(s.D);
                            d2=testschnitzcells(s.E);
                            mutimeidxd1=find(d1.(muTimeField)>MuTime-timeTolerance &  d1.(muTimeField)<MuTime+timeTolerance);
                            mutimeidxd2=find(d2.(muTimeField)>MuTime-timeTolerance &  d2.(muTimeField)<MuTime+timeTolerance); 
                            % if both idx are not empty (daughters live long enough)
                            if ~isempty(mutimeidxd1) & ~isempty(mutimeidxd2) 
                               mu1=d1.muP15_fitNew_cycCor(mutimeidxd1)*d1.birthLengthRatio15+ ... 
                                        d2.muP15_fitNew_cycCor(mutimeidxd2)*d2.birthLengthRatio15; %average growth rates
                               if ~isnan(p1+mu1)
                                  mup_vec1=[mup_vec1; mu1, p1,i];
                               end    
                            else
                              %  disp(['Daughter cells (at least one) don''t contain data at right time distance for schnitz ' num2str(i)]);
                            end
                        end
                % ---------------------------------------
                % Elseif: Schnitz has no data of Mu. Mu is earlier, search parent
                % ---------------------------------------        
                    else
                        if s.P~=0
                            par=testschnitzcells(s.P);
                            mutimeidxpar=find(par.(muTimeField)>MuTime-timeTolerance &  par.(muTimeField)<MuTime+timeTolerance);
                            % if idx is not empty (parent got old enough)
                            if ~isempty(mutimeidxpar)
                               mu1=par.muP15_fitNew_cycCor(mutimeidxpar); % use growth rate directly
                               if ~isnan(p1+mu1)
                                  mup_vec1=[mup_vec1; mu1, p1,i];
                               end    
                            else
                               % disp(['Parent cell doesn''t contain data at right time distance for schnitz' num2str(i)]);
                            end
                        end
                    end
                end
                % --- end of: loop over different cases whether muTims is at
                % schnitz, P or D&E
            end
            % --- end of: time points within fitTime?
        end
        % --- end of: loop over time points
    end
    % --- end: if useForPlot
end
% end of: loop over schnitzes
                    
disp(['Number of datapoints: ' num2str(size(mup_vec1,1))]);            
            
lastone=length(mup_vec1);

x=mup_vec1(1:lastone,1);
y=mup_vec1(1:lastone,2);
schn=mup_vec1(1:lastone,3);
clf
plot(x,y,'.')
                
corrcoef(x,y)                
            
        
    
    