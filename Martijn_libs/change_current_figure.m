function change_current_figure(h)
% Change figure without getting it into focus
% Thanks to:
% http://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
set(0,'CurrentFigure',h)