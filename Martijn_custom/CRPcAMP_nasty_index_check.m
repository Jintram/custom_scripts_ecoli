
FIELD1='time_MW_atdX';
FIELD2='dG5_cycCor';
%FIELD2='muP15_fitNew_all_MW_atdX';


disp(['***' 10 'Starting analysis.']);

numberofLengthIssues=0; numberofValueIssues=0;
for i = 1:numel(schnitzcells)
    
    if (numel(schnitzcells(i).(FIELD1)) ~= numel(schnitzcells(i).(FIELD2)))
        deltaNumel = numel(schnitzcells(i).(FIELD1)) - numel(schnitzcells(i).(FIELD2));
        disp(['There are issues with length of fields being consistent at schnitz ' num2str(i) ' (delta=' num2str(deltaNumel) ')']);
        numberofLengthIssues=numberofLengthIssues+1;
    elseif ~any(size(schnitzcells(i).(FIELD1))==0)
        if(any( schnitzcells(i).(FIELD1) ~= schnitzcells(i).(FIELD2) ))
            %disp('There are issues with values being consistent.');
            numberofValueIssues=numberofValueIssues+1;
        end
    end
        
end

disp(['***' 10 'There are ' sprintf('%5.f',numberofLengthIssues) ' length issues..' ...
        10 'There are ' sprintf('%5.f',numberofValueIssues) ' value issues..' ]);
