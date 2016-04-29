% Adds the frame number to each schnitzcells in which a certain
% fluorescence production rate (e.g. dY5) exists. Similar to 'Y_frames' for
% Y-concentration. It's necessary for applying NW_plot_cellCycleDependence
% and cell cycle corrections (these functions rely on frame field and not
% on time field)
% Necessary input: timefield at which the fluor data exists
% This input field will be compared to 'time' (existent for each phase
% image)
% The frame numbers are then copied from 'frames' for the appropriate time
% points
%
%
% ********* ADJUST **********
myschnitzcells=schnitzcells;   % FUNCTION DOES NOT TEST FOR USE FOR PLOT!!!  AND SCHNTIZCELLS MUST BE REASSOCIATED AT END OF FCT!!!
fluotime='dC5_time';  % name of input field
framefield='dC5_frames';    % name of newly made frame field

timetolerance=0.5;  % since comparing 'double' values, allow for some tolerance and don't inforce ==

for i=1:length(myschnitzcells)
    % ************* NO USEFORPLOT CONTROL!!! **********
    myschnitzcells(i).(framefield)=[];   % initiate new field
    s=myschnitzcells(i);
    if ~isempty(s.(fluotime))   % fluor rate data for at leat one time point
        for ff=1:length(s.(fluotime))
            idx=find(s.time>s.(fluotime)(ff)-timetolerance & s.time<s.(fluotime)(ff)+timetolerance);
            if length(idx)~=1
                disp(['Don''t find appropriate time data for schnitz ' num2str(i) ', or found too many.']);
                break
            else
                s.(framefield)(end+1)=s.frames(idx);
            end
        end
    end
    myschnitzcells(i)=s;
end




 schnitzcells=myschnitzcells;