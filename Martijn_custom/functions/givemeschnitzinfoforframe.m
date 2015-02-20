
function [theSchnitzes,cellTimePoints,cellLengths,timemin,timemax,lengthmin,lengthmax] =  givemeschnitzinfoforframe(p, schnitzcells, fr)
% function [cellTimePoints,cellLengths] =  givemelengthsforframe(p, schnitzcells, fr)
% Returns [time, lengths] of all schnitzes in the given frame, at the
% time of that frame.

cellTimePoints=[]; cellLengths=[]; theSchnitzes=[];
for i = 1:numel(schnitzcells)
    
    % backwards compatibility
    if isfield(schnitzcells,'frames')
        schnitzcells(i).frame_nrs = schnitzcells(i).frames-1; % N+1 bug
    end
    
    % find if and where this schnitzcell contains info about the current
    % frame
    hit = find(fr==schnitzcells(i).frame_nrs);
    if ~isempty(hit)        
                
        % store that data in an array
        cellTimePoints(end+1) = schnitzcells(i).time(hit);
        cellLengths(end+1)    = schnitzcells(i).length_fitNew(hit); % check whether this is correct length
        theSchnitzes(end+1)   = i;
    end
    
    % Note that this code is redundant:
    timemin=min([schnitzcells(:).time]);
    timemax=max([schnitzcells(:).time]);
    lengthmin=min([schnitzcells(:).length_fitNew]);
    lengthmax=max([schnitzcells(:).length_fitNew]);
    
    
end

% plot data
%{
figure, hold on
datax={[1,1],[2,2]}
datay={[2,3],[3,6]}
for i =1:numel(datax)
    plot(datax{i}, datay{i},'o')
end
xlim([0,10])
ylim([0,10])
%}


