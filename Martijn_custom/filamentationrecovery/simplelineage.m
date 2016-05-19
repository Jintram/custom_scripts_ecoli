


endSchnitz = numel(schnitzcells);

currentSchnitz=endSchnitz;
reverseLineageList=[]; count=0;
while 1
   
    parentSchnitz = schnitzcells(currentSchnitz).P;
    
    if parentSchnitz == 0
        break;
    end
    
    reverseLineageList(end+1) = parentSchnitz;
    currentSchnitz = parentSchnitz;
    
    % infinite loop protection
    count = count+1;
    if count>10000, error('Broke while loop to avoid infinite loop.'); end
end

lineageList = reverseLineageList(end:-1:1);