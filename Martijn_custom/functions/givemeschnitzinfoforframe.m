
function [theSchnitzes,cellTimePoints,cellLengths,timemin,timemax,lengthmin,lengthmax,mapIndexes,output] =  givemeschnitzinfoforframe(p, schnitzcells, fr)
% function [theSchnitzes,cellTimePoints,cellLengths,timemin,timemax,lengthmin,lengthmax] =  givemeschnitzinfoforframe(p, schnitzcells, fr)
%
% Returns [time, lengths] of all schnitzes in the given frame, at the
% time of that frame.
% INPUT:
% p                 our standard input parameter
% schnitzcells      our standard schnitz datastruct
% fr                framenumber you're interested in
%
% OUTPUT
% theSchnitzes      schnitzes that are member of this frame
% cellTimePoints    the timepoints at which the different found schnitzes live
% cellLengths       idem for lenghts
% timemin           min value of cellTimePoints
% timemax           max ""
% lengthmin         min value of cellLengths
% lengthmax         max ""
% mapIndexes        the indexes of the schnitzes in the current frame

cellTimePoints=[]; cellLengths=[]; theSchnitzes=[]; mapIndexes=[];
output = {}; totalHits = 0;
for i = 1:numel(schnitzcells)
    
    % backwards compatibility
    if isfield(schnitzcells,'frames')
        schnitzcells(i).frame_nrs = schnitzcells(i).frames-1; % N+1 bug
    end
    
    % find if and where this schnitzcell contains info about the current
    % frame
    hit = find(fr==schnitzcells(i).frame_nrs);
    if ~isempty(hit)        
                
        totalHits = totalHits + 1;
        
        % store that data in an array
        cellTimePoints(end+1) = schnitzcells(i).time(hit);
        cellLengths(end+1)    = schnitzcells(i).length_fitNew(hit); % check whether this is correct length
        theSchnitzes(end+1)   = i;
        mapIndexes(end+1)     = schnitzcells(i).cellno(hit);
        
        % all info that's possible for that schnitz:
        myFieldnames = fieldnames(schnitzcells(i));
        for fieldno=1:numel(myFieldnames)
            if ismatrix(schnitzcells(i).(myFieldnames{fieldno})) && numel(schnitzcells(i).(myFieldnames{fieldno}))>=hit % note this is somewhat pathological when numel = 1            
                output{totalHits}.(myFieldnames{fieldno}) = schnitzcells(i).(myFieldnames{fieldno})(hit);
            end
        end
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


