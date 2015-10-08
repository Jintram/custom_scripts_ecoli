function count = MW_countcellsinlastframe(schnitzcells)
% A super simple function that identifies the last frame and then brute
% force counts the number of schnitzcells that are associated to that
% frame.

count=0;

frames = unique([schnitzcells.frame_nrs]);
lastFrame=frames(end);

for i = 1:numel(schnitzcells)
    if any(schnitzcells(i).frame_nrs == lastFrame)
        count = count+1;
    end
end

%count

end