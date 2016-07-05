
%% Checks for inconsistencies in length of each schnitz infos.

for i = 1:numel(schnitzcells)
    
    if (numel(schnitzcells(i).frame_nrs))~=(numel(schnitzcells(i).ceny))
        disp(num2str(i)); 
    end
    
end

%% check for NaN values

FIELDOFINTEREST = 'G6_mean';
FIELDOFINTEREST = 'schnitzcells.muP5_fitNew_atG6';
FIELDOFINTEREST = 'muP5_fitNew_cycCor';


disp('NaN values:');

for i = 1:numel(schnitzcells)
    if any( isnan(schnitzcells(i).(FIELDOFINTEREST)) );
        disp(num2str(i)); 
    end
end