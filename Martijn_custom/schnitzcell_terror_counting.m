
% Checks for inconsistencies in length of each schnitz infos.

for i = 1:numel(schnitzcells)
    
    if (numel(schnitzcells(i).frame_nrs))~=(numel(schnitzcells(i).ceny))
        disp(num2str(i)); 
    end
    
end