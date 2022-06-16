function [schnitz2]=NW_addToSchnitzes_FluorIndex(p,colorfield)
% This function adds an extra Field for each schnitz and each time point
% with fluorescence data. The field is an index telling on which position
% of all fluorescent data at each time point the schnitz lies.
% At each time point, the brightest Schnitz has index 1, the darkest
% Schnitz index '#schnitzes at this time point'. Is in the spirit of
% Rosenfeld2005
%
% MAYBE RENUMBERING OF INDEX IS NECESSARY WHEN SOME SCHNITZES ARE REMOVED
% (or check useForPlot... and use a 'schntizcell' as fct input and don't
% load it)
%
% vielleicht sollte der Index auch erst bei den Branches hinzugefuegt
% werden!!
%
% FUNCTION IS NOT YET FOOLPROOF!

% load schnitzcells (only default folder possible)
schnitzname = [p.tracksDir,p.movieName,'-Schnitz.mat'];
load(schnitzname);

% create new variable name
fluorindex=[upper(colorfield) '_index'];
fluorfield=[upper(colorfield) '6_mean'];

% initiate index fields with Nan (assures at least that field has same
% length as fluorfield
for i=1:length(schnitzcells)
    schnitzcells(i).(fluorindex)=NaN(1,length(schnitzcells(i).(fluorfield)));
end    

% get timeField and unique time points
timefield = [upper(colorfield) '_time'];
allTimePoints=[];
for i=1:length(schnitzcells)
    allTimePoints=[allTimePoints, schnitzcells(i).(timefield)];
end
unique_timeField  = unique(allTimePoints);

% loop over all time points
for tt=unique_timeField
    schnitzFluoArray=[];
    % loop over all schnitzes
    for schnitz=1:length(schnitzcells) % useForPlot could be here included
        %check if schnitz exist at time tt and if yes record schnitz number
        %and fluorescence
        idx=find(schnitzcells(schnitz).(timefield)==tt);
        if ~isempty(idx)
            schnitzFluoArray=[schnitzFluoArray; schnitz, schnitzcells(schnitz).(fluorfield)(idx)];
        end   
    end
    
    % sort array in descending order
    arraySorted=sortrows(schnitzFluoArray,-2);  %-2: 2nd column, descending
    
    %if tt>100 & tt<250  %debugging
    %    tt
    %    arraySorted
    %end 
    
    % write index into schnitz info
    for i=1:size(arraySorted,1)
        schn=arraySorted(i,1);
        idx2=find(schnitzcells(schn).(timefield)==tt); %must be non-empty
        schnitzcells(schn).(fluorindex)(idx2)=i;
    end
end

schnitz2=schnitzcells;
    save(schnitzname,'schnitzcells');