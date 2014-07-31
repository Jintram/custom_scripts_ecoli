% adds logarithm to length field

% 

%log field: 2-log of length
%logdivPhasefield: first value of log field substracted from logFieldArray,
%   then divided by phase -> if exponential growth, it should result in a
%   line: length=length_init*2^(time/dbltime);
%logdivTtimefield: time calculated as in NW_plotCellcycle dependence
%logdivTinvtimefield: time calculated as in NW_plotCellcycle dependence.
%reference length is last measured length (contrary to first in case above)

%****ADJUST*****
%lengthfield='av_fittedLength_fitNew';
lengthfield='length_fitNew';
myschnitzcells=s_rm_fitTime_cycle;
%***************

% new fields
logfield=['log_' lengthfield];
logdivPhasefield=['logdivPhase_' lengthfield];
logdivPhaseinvfield=['logdivPhaseinv_' lengthfield];

logdivTfield=['logdivT_' lengthfield];
logdivTinvfield=['logdivTinv_' lengthfield];

for i=1:length(myschnitzcells)
    s=myschnitzcells(i);
    
    s.(logfield)=log(s.(lengthfield))./log(2);
    
    s.(logdivPhasefield)=s.(logfield)-s.(logfield)(1);
    s.(logdivPhasefield)=s.(logdivPhasefield)./s.phase;
    
    s.(logdivPhaseinvfield)=s.(logfield)-s.(logfield)(end);
    s.(logdivPhaseinvfield)=s.(logdivPhaseinvfield)./(s.phase-s.phase(end));
    
    totalTime=length(s.frames); %# frames a schnitz lives
    YdataTime=s.frames-s.frames(1); %absolute time with initial frame set to 0
    YdataTime=YdataTime/totalTime; % relative time (birth=0. last frame=(#frames-1)/#frames
    YdataTimeInv=YdataTime-YdataTime(end);
    
    s.(logdivTfield)=s.(logfield)-s.(logfield)(1);
    s.(logdivTfield)=s.(logdivTfield)./YdataTime;
    
    s.(logdivTinvfield)=s.(logfield)-s.(logfield)(end);
    s.(logdivTinvfield)=s.(logdivTinvfield)./YdataTimeInv;
    
    myschnitzcells(i).(logfield)=s.(logfield);
    myschnitzcells(i).(logdivTfield)=s.(logdivTfield);
    myschnitzcells(i).(logdivPhasefield)=s.(logdivPhasefield);
    myschnitzcells(i).(logdivTinvfield)=s.(logdivTinvfield);
     myschnitzcells(i).(logdivPhaseinvfield)=s.(logdivPhaseinvfield);
    
    %ff=['bla'];
    %myschnitzcells(i).(ff)=s.(logfield)-s.(logfield)(end);
    %disp([num2str(myschnitzcells(i).(ff)(end))]);
end

schnitzcells_length=myschnitzcells;

clear myschnitzcells s logfield logdivTfield YdataTime totalTime YdataTimeInv;

    