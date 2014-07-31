% creates a histogram of the # phase contrast images (frames) that are
% typically obtained within the cell cycle of one bacteria. Can be used to
% determine fitting window of growth rate.

% if useForPlot is set, only cells which have =1 will be considered.

myschnitzcells=s_rm_fitTime;
myschnitzcellsname='s_rm_fitTime';



framesvec=[];
if isfield(myschnitzcells,'useForPlot')
    PLOTSET=1;
else
    PLOTSET=0;
end
    
for i=1:length(myschnitzcells)
    if PLOTSET
        if (myschnitzcells(i).useForPlot==1 & myschnitzcells(i).completeCycle==1)
            framesvec=[framesvec; length(myschnitzcells(i).frames)];
        end
    else
        if (myschnitzcells(i).completeCycle==1)
            framesvec=[framesvec; length(myschnitzcells(i).frames)];
        end
    end
end

figure
clf
hist(framesvec)
title(['# frames per cecy: mean=' num2str(mean(framesvec)) ', median=' num2str(median(framesvec)) '. Used ' myschnitzcellsname],'Interpreter','None')
xlabel('# frames per cecy')
ylabel('frequency')
    