% calculates the average number of frames that a cell exists -> helpful to
% determine which averaging frame size XX in muPXX_fitNew should be used
% (rule of thumb: 1/3 of cell cycle)

myschnitzcells=s_rm_fitTime; % use one selection where 'weird' slow schnitzes etc are removed and mu is cst

frames_unique=[];
for i=1:length(myschnitzcells)
    s=myschnitzcells(i);
    if s.useForPlot==1 & s.completeCycle==1     
       frames_unique=[frames_unique; length(s.frames)];
    end
end

meanframes=mean(frames_unique);
stdframes=std(frames_unique);

disp('  ')
disp('-----------------------')
disp(['Cells exist on average ' num2str(meanframes) ' frames.  Stdev is ' num2str(stdframes) '.'])
disp('-----------------------')