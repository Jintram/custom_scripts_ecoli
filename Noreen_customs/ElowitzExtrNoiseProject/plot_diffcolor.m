function []=plot_diffcolor(vecX,vecY,matColor,vecMarker)
%plots vecY vs vecX with each datapoint in a different color, specified by
%matColor and symbol kind of specified by vecMarker
%
% figure should be opened & clf beforehand
%
% ------
% Input
% ------
% vecX: nx1 vector, x-axis
% vecY: nx1 vector: y-axis
% matColor: nx3 matrix. Color matrix: each entry btw 0&1
% vecMarker: nx1 vector (cell array) specifying the symbol (only distinguishes
%                        btw '.' and 's' and other)

if length(vecX)~=length(vecY) | length(vecX)~=size(matColor,1) ... 
        | length(vecX)~=length(vecMarker)
    error('Error: vectors & color matrix must be same length.')
end

hold on
for i=1:length(vecX)
    if strcmp(vecMarker(i),'.')==1
        plot(vecX(i),vecY(i),'.','MarkerSize',25,'Color',matColor(i,:))
    elseif strcmp(vecMarker(i),'s')==1
        plot(vecX(i),vecY(i),'s','MarkerSize',8,'Color',matColor(i,:),'MarkerFaceColor',matColor(i,:))    
    else
        plot(vecX(i),vecY(i),'^','MarkerSize',8,'Color',matColor(i,:),'MarkerFaceColor',matColor(i,:))    
    end
end