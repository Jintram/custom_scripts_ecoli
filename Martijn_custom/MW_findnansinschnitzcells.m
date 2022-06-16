
%% Can tell you whether fields contain NaN values in schnitzcells structure

%% 

FIELDOFINTEREST = 'muP5_fitNew_all';

for i = 1:numel(schnitzcells)
    
    if any(isnan(schnitzcells(i).(FIELDOFINTEREST)))
        disp([num2str(i) ' is suspicious']);
    end
    
end