function [myschnitzcells] = NW_normalizeFluor_addDifference(p,myschnitzcells,Field1,Field2,Timefield,DiffFieldName) 
% This function normalizes each fluorescent value (field1, field2) by the colony average at
% this time point. (Alternately the total over-time average could also be
% taken but I (NW) thought that global temporal trends might make this
% latter optinon unsuitable). The normalized fields are called Field1_N and
% Field2_N.
% Additionally, a field with the difference of the normalized fields (i.e.
% Field1_N - Field2_N)  is added. This fields represents intrinsic noise. 
%
% Note that I find it more suitable to only take healthy cells into account
% for calculating the average (since sick cells often have abnormally high
% or low fluorescence). Therefore, a "myschnitzcells" should be given where
% useForPlot is set (e.g. s_rm)
% For test reasons, also fields that are normalized with the total average
% over all time points and all (useForPlot==1) cells, are added:
% field1_Ntot, field2_Ntot
%
% ********************************************************************
% THERE IS NO (NOT YET) FANCY CHECK FOR EVERY POSSIBLE INPUT ERROR!
% SO GIVE THE CORRECT INPUTS!
% ********************************************************************
% PAY ATTENTION IF ABSOLUTE (Betrag) or DIFFERENCE NOISE IS ADDED!!!
%
% CURRENTLY IGNORES MYSCHNITZCELLS INPUT AND LOADS SCHNITZCELLS_RM FILE!!!!! (MUST
% ALREADY EXIST ....very akward
% *******************************************************************
%
% REQUIRED INPUT
% - Field1, Field2: fluorescence fields. e.g. Y6_mean, C5_mean,
% dY5_sum_dt_subtr. They have to have equal length if you don't want the
% program to get lost.
% - myschnitzcells: best use s_rm (or s_rm_fitTime)
% - p
% - DiffFieldName: how the difference field (field1_N-field2_N) is supposed
%             to be named
% - Timefield: timefield corresponding to Field1/Field2. Typically: Y_time
%
% OUTPUT
% - myschnitzcells: schnitzcells with added fields. Is also being saved in
%                  posXcrop/data/ folder in posXcrop-Schnitz-DiffAdded.mat


% **************IMPROVE!!!****************************
%schnitzname = [p.tracksDir,p.movieName,'-Schnitz-DiffAdded.mat']; %blubb
%load(schnitzname); %blubb
%myschnitzcells=schnitzcells_rm; %blubb
% ***********************************************

% new field names
Field1Norm=[Field1, '_N']; % divided by each time point mean
Field2Norm=[Field2, '_N'];
Field1NormTot=[Field1, '_Ntot'];
Field2NormTot=[Field2, '_Ntot'];
DiffFieldNametot=[DiffFieldName, 'tot'];


% get unique array of time points
timearray=[];
for i=1:length(myschnitzcells)
    if length(myschnitzcells(i).(Timefield))>0
        timearray=[timearray, myschnitzcells(i).(Timefield)];
    end
end
timearray=unique(timearray);


% ----------------------------------------------------------
% get average value for each time point
% ----------------------------------------------------------
%
meanField1=['mean' , Field1];
meanField2=['mean' , Field2];
%
data=struct;
% data(k) corresponds to timepoint "timearray(k)". it contains the actual
% time and an array of all Field1 and Field2 values existing at this time
% point.
% necessary? -> preallocation
for i=1:length(timearray);
    data(i).time=timearray(i);
    data(i).(Field1)=[];
    data(i).(Field2)=[];
    data(i).(meanField1)=[];
    data(i).(meanField2)=[];
end

%loop over schnitzes
for run=1:length(myschnitzcells)
    
    s=myschnitzcells(run);
    if s.useForPlot==1
        %loop over each time points in which schnitz exists and add to
        %"data"
        if length(s.(Field1))>0 & length(s.(Timefield))>0 %dummy check for stupid data
            for f=1:length(s.(Field1));
                mytime=s.(Timefield)(f); % should always work since time field is same length or 1 longer than Field1
                idx=find(timearray==mytime);
                if ~isnan(s.(Field1)(f)) &  ~isnan(s.(Field2)(f)) %should always be true but you never know...
                    data(idx).(Field1)=[data(idx).(Field1), s.(Field1)(f)];
                    data(idx).(Field2)=[data(idx).(Field2), s.(Field2)(f)];
                end
            end
        end
    end
end
% total mean (cells equally weighed)
field1all=[];
field2all=[];
%add extra fields to data with average for each time points
for i=1:length(data)
    data(i).(meanField1)=mean(data(i).(Field1));
    data(i).(meanField2)=mean(data(i).(Field2));
    field1all=[field1all,data(i).(Field1)];
    field2all=[field2all,data(i).(Field2)];
    
end
meanOfEverythingField1=mean(field1all);
meanOfEverythingField2=mean(field2all);

% add normalized fields and difference field to schnitzcells
for i=1:length(myschnitzcells)
    % assign empty arrays as a start
    myschnitzcells(i).(Field1Norm)=[];
    myschnitzcells(i).(Field2Norm)=[];
    myschnitzcells(i).(Field1NormTot)=[];
    myschnitzcells(i).(Field2NormTot)=[];
    myschnitzcells(i).(DiffFieldName)=[];
    myschnitzcells(i).(DiffFieldNametot)=[];
    clear s;
    s=myschnitzcells(i);
    % assign real data
    if s.useForPlot==1
        if length(s.(Field1))>0 & length(s.(Timefield))>0 %dummy check for stupid data
            for f=1:length(s.(Field1))
                mytime=s.(Timefield)(f); % should always work since time field is same length or 1 longer than Field1
                idx=find(timearray==mytime);
                s.(Field1Norm)(f)=s.(Field1)(f)/data(idx).(meanField1);
                s.(Field2Norm)(f)=s.(Field2)(f)/data(idx).(meanField2);
                s.(Field1NormTot)(f)=s.(Field1)(f)/meanOfEverythingField1;
                s.(Field2NormTot)(f)=s.(Field2)(f)/meanOfEverythingField2;
                s.(DiffFieldName)=s.(Field1Norm)-s.(Field2Norm); %blubb
                s.(DiffFieldNametot)=s.(Field1NormTot)-s.(Field2NormTot);%blubb
                %s.(DiffFieldName)=abs(s.(Field1Norm)-s.(Field2Norm));%blubb
                %s.(DiffFieldNametot)=abs(s.(Field1NormTot)-s.(Field2NormTot));%blubb
            end
        end
    end
    myschnitzcells(i)=s;
end

%save schnitzcells under new name
%schnitzcells_rm=myschnitzcells;        %blubb
%schnitzname = [p.tracksDir,p.movieName,'-Schnitz-DiffAdded.mat']; %blubb
%save(schnitzname,'schnitzcells_rm'); %blubb
