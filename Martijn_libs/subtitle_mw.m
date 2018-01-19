function [ax,h]=subtitle_mw(text,theXLabel,theYLabel,locx1,locy1)
% function [ax,h]=subtitle_mw(text,theXLabel,theYLabel)
%
% SOURCE: https://nl.mathworks.com/matlabcentral/answers/100459-how-can-i-insert-a-title-over-a-group-of-subplots
% USER:   MathWorks Support Team
%
%Centers a title over a group of subplots.
%
%Returns a handle to the title and the handle to an axis.
%
% [ax,h]=subtitle(text)
%
%           returns handles to both the axis and the title.
%
% ax=subtitle(text)
%
%           returns a handle to the axis only.

if ~exist('locx1','var')
    locx1=.075;
end
if ~exist('locy1','var')
    locy1=.05;
end

% Create an invisible axis that covers all the subplots
ax=axes('Units','Normal','Position',[locx1 locy1 .85 .9],'Visible','off');
%ax=axes('Units','Normal','Position',[.075 .075 .85 .85]);%,'Visible','off'); % original

% Remove ticklabels
set(ax, 'XTickLabel', []);
set(ax, 'YTickLabel', []);

% Set desired components to visible
set(get(ax,'Title'),'Visible','on')
set(get(ax,'xlabel'),'Visible','on')
set(get(ax,'ylabel'),'Visible','on')

% print the desired components
title(text);
xlabel(theXLabel);
ylabel(theYLabel);

% some additional output
if (nargout < 2)

    return

end

h=get(ax,'Title');
