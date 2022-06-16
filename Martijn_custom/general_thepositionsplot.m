
% obtain positions
% ===

% position file
%myFile = 'F:\A_Tans0_step1_incoming_not_backed_up\2014-01-31\positionlist_2014-01-31.pos';
myFile = 'F:\A_Tans0_step1_incoming_not_backed_up\2014-01-23 - evening\positionlist_2014-01-23_evening.pos';
%myFile = 'F:\A_Tans2_step4c_failed_entirely\2014-01-24 - morning\positionlist_2014-01-24_morning.pos'

% read position file
contents = fileread(myFile);
contents = strsplit(contents,'\n'); % process linebreaks

xvalues=[];
yvalues=[];

% loop over file contents per line
for i=1:numel(contents)

    % if line for stage position found
    if strcmp(contents{i},'               "DEVICE": "XYStage",')    
        
        % load x and y lines (located at +3 and +3 resp. for X,Y)
        xstr = contents{i+3};
        ystr =contents{i+2};
        
        % and extract position values for x and y, save to list
        expression = '[ -]\d*\.\d*|[ -]\d+';
        xvalues(end+1) = str2num(cell2mat(regexp(xstr,expression,'match')));
        yvalues(end+1) = str2num(cell2mat(regexp(ystr,expression,'match')));
        
    end
    
end

% make plot
% ===          

figure, plot(xvalues,yvalues,'o','MarkerSize',7,'LineWidth',3);
xlabel('X');
ylabel('Y');
title(['Pad positions (folder numbering used)', 10, myFile]);

% labels
for i=[1:numel(xvalues)]
    text(xvalues(i),yvalues(i),['pos' num2str(i)],'Color',[.5 .5 .5]);
end

% print for convenience file name
myFile



