% extracts area, total fluorescence and fluorescence rates and
% concentrations frrom schnitzcells
% calculates relative increase for each time step

% **** ADJUST *****
schnitzUse=s_rm_fitTime;  %schnitzcells selection
yfprate='dY5_sum_dt_subtr';
cfprate='dC5_sum_dt_subtr';
area='area';
yfptot='Y5_sum';
cfptot='C5_sum';
yfpconc='Y5_mean';
cfpconc='C5_mean';
% *****************

Area=[];
YFPrate=[];
CFPrate=[];
YFPtot=[];
CFPtot=[];
YFPconc=[];
CFPconc=[];


% loop over all schnitzes
for i=1:length(schnitzUse)
    s=schnitzUse(i);
    % use only schnitzes with useForPlot=1
    if s.useForPlot==1
        % check if fields actually contain data (length>0)
        % Note: e.g. prod.rates and mu do not always have same array length
        % (rates can be 1 shorter for late data)
        minLength=min([length(s.(field1)),length(s.(field2)),length(s.(field3))]);
        if minLength>0
            % check if NaN values exist
            Nan1=sum(isnan(s.(field1)(1:minLength)));
            Nan2=sum(isnan(s.(field2)(1:minLength)));
            Nan3=sum(isnan(s.(field3)(1:minLength)));
            if (Nan1+Nan2+Nan3)==0
                % add values to mymatrix
                addmatrix=[s.(field1)(1:minLength)',s.(field2)(1:minLength)',s.(field3)(1:minLength)'];
                mymatrix=[mymatrix;addmatrix];
            end
        end
    end
end