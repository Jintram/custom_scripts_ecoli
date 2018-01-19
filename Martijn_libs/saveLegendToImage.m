function saveLegendToImage(figHandle, legHandle, subHandles, fileName, fileType)
% source: https://stackoverflow.com/questions/18117664/how-can-i-show-only-the-legend-in-matlab
% written by user Franz Wurst
%
% Edits to handle subplots and not write image by MW
% Set subHandles to empty matrix [] if it's a figure w/o subplots

%make all contents in figure invisible
allLineHandles = findall(figHandle, 'type', 'line');

for i = 1:length(allLineHandles)

    allLineHandles(i).XData = NaN; %ignore warnings

end

%make axes invisible
axis off
% In case there are subplots (addition by MW)
for hIdx=1:numel(subHandles)
    subplot(subHandles(hIdx));    
    axis off;
    title([]);
end

%move legend to lower left corner of figure window
legHandle.Units = 'pixels';
boxLineWidth = legHandle.LineWidth;
%save isn't accurate and would swallow part of the box without factors
legHandle.Position = [6 * boxLineWidth, 6 * boxLineWidth, ...
    legHandle.Position(3), legHandle.Position(4)];
legLocPixels = legHandle.Position;

%make figure window fit legend
figHandle.Units = 'pixels';
figHandle.InnerPosition = [1, 1, legLocPixels(3) + 12 * boxLineWidth, ...
    legLocPixels(4) + 12 * boxLineWidth];

%save legend
if exist('fileName','var') % if statement = addition MW to make it optional
    saveas(figHandle, [fileName, '.', fileType], fileType);
end

end