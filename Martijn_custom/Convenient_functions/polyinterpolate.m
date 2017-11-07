function [fittedLineX,fittedLineY] = polyinterpolate(xval,yval,windowWidth,pointsPerSection)

% fit settings
% windowWidth=3; 
% pointsPerSection = 20;

% distance between two fitted points, based equidistant input
dx=(xval(2)-xval(1))/pointsPerSection; % points should be equidistant
% complete set of x values for extrapolation
fittedLineX = min(xval):dx:max(xval); 
% set up y-line
fittedLineY= zeros(1,numel(fittedLineX));
% set a weighing such that sections are smoothly merged
weighing        = [0:pointsPerSection pointsPerSection-1:-1:0]/pointsPerSection;
weighingstart   = [zeros(1,pointsPerSection) pointsPerSection:-1:0]/pointsPerSection;
weighingend     = [0:pointsPerSection zeros(1,pointsPerSection)]/pointsPerSection;

for windowIdx=windowWidth+1:numel(xval)-windowWidth
    % get current points to fit to 
    currentxvaltofit=xval(windowIdx-windowWidth:windowIdx+windowWidth);
    currentyvaltofit=yval(windowIdx-windowWidth:windowIdx+windowWidth);

    % perform the fit
    p=polyfit(currentxvaltofit,currentyvaltofit,2);

    % calculated the line from fit
    xvalfitline = xval(windowIdx-1):dx:xval(windowIdx+1);
    yvalfitline = p(1).*xvalfitline.^2 + p(2).*xvalfitline + p(3);

    % 
    sectionIndices=(1+pointsPerSection*(windowIdx-2)):(1+pointsPerSection*(windowIdx));
    if windowIdx==windowWidth+1 
        % special case 1st point
        fittedLineY(sectionIndices) = fittedLineY(sectionIndices)+yvalfitline.*weighingstart;
    elseif windowIdx==numel(xval)-windowWidth
        % special case last point
        fittedLineY(sectionIndices) = fittedLineY(sectionIndices)+yvalfitline.*weighingend;
    else            
        % every other point
        fittedLineY(sectionIndices) = fittedLineY(sectionIndices)+yvalfitline.*weighing;
    end

end

determined = (fittedLineY>0);

fittedLineX=fittedLineX(determined);
fittedLineY=fittedLineY(determined);


end