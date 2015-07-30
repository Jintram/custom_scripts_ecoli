
function slide_plotmore(p,schnitzcells,ts)
% Using Daan's function DJK_analyzeMu, plots with fits are produced 
% automatically for different fit windows.
%
% provide p, schnitzells as input arguments, and also ts, which can be
% an array with all timepoints, but also just contain the maximum time.
%
% further configuration: nrtbins can be set in function

    nrtbins=10
    maxt=max(ts);
    mytwindows=[0:maxt/nrtbins:maxt];

    % whole range
    fitTime = DJK_analyzeMu(p, schnitzcells, 'xlim', [0 maxt], 'onScreen', 0,'fitTime',[0 maxt]);

    % subranges
    for idx_window = 1:(length(mytwindows)-1)

        fitTime = DJK_analyzeMu(p, schnitzcells, 'xlim', [0 maxt], 'onScreen', 0,'fitTime',[mytwindows(idx_window) mytwindows(idx_window+1)]);

    end

end