
%% Checks for inconsistencies in length of each schnitz infos.

for i = 1:numel(schnitzcells)
    
    if (numel(schnitzcells(i).frame_nrs))~=(numel(schnitzcells(i).ceny))
        disp(num2str(i)); 
    end
    
end

%% check for NaN values
FIELDSOFINTEREST = {'G6_mean','muP5_fitNew_cycCor','dG5_cycCor'};
FIELDSOFINTEREST = {'C6_mean','muP5_fitNew_cycCor','dC5_cycCor','Y6_mean','dY5_cycCor'};

disp('===');
listAllBadOnes=[];
for idx=1:numel(FIELDSOFINTEREST)
    
    FIELDOFINTEREST=FIELDSOFINTEREST{idx};

    disp(['NaN values for ' FIELDOFINTEREST ':']);

    listOfBadOnes=[];
    for i = 1:numel(schnitzcells)
        if any( isnan(schnitzcells(i).(FIELDOFINTEREST)) );
            %disp(num2str(i)); 
            listOfBadOnes(end+1)=i;
        end
    end

    disp(num2str(listOfBadOnes));
    listAllBadOnes=[listAllBadOnes listOfBadOnes];
    
end
listAllBadOnes=unique(listAllBadOnes);
disp(['All: ' 10 num2str(listAllBadOnes)]);



