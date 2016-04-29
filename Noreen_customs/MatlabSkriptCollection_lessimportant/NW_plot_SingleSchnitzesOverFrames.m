function [alldata,myfieldmeantrace]=NW_plot_SingleSchnitzesOverFrames(p,myschnitzcells,myfield,myframes)

% plots traces of single schnitzes over time and their average.
% IT'S NOT FOOL PROOF!
% currently, 'frames' or 'R_frames'  is the only option for x-axis input. i.e. frames is
% the parameter for time.
% NO WEIGHING WHATSOEVER IS PERFORMED. ONLY DIVIDED BY THE NUMBER OF
% SCHNITZES FOR AVERAGE
% THIS FUNCTION IS MAINLY INTENDED FOR BEADS THAT EXIST OVER THE WHOLE
% MOVIE TIME.
%
% Required arguments:
% p:             stuct. needed to determine save-folder
% myschnitzcells:  structure fom 'fullAnalysis' which has complete tracking
%                data
% myfield:       which field to plot at y-axis (e.g. R6_mean, area)
% myframes:      "time" axis, use either 'frames' or 'R_frames'



% to come: chekc for all possible stupid inputs one can give

% to come: allow for other input-x-axis than 'frames'


%make save directory
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'schnitzcells' filesep];
end
% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end

%--------------------------------------------------------------------------
% Collect all relevant data and plot
%--------------------------------------------------------------------------

alldata=zeros(0,3);
% contains complete x&y data of all schnitzes
% schnitznum1  ---  x (framenr1) --- myfield1
% schnitznum1  ---  x (framenr2) --- myfield2
% ...
% schnitznum2  ---  x (framenr2) --- myfield2
% schnitznum2  ---  x (framenr2) --- myfield2
% ...

% e.g.
% 1 --- 8 (first frame nr) --- 444 (conc)
% 1 --- 16 (seconde frame nr) --- 445 (conc)
% ...
% 1 --- 80 (last frame nr) --- 440 (conc)
% 2 --- 16 (first frame nr) --- 821 (conc)
% ...

fig1=figure(20);
set(gcf,'WindowStyle','docked')
clf
hold on
xlabel('frame number')
ylabel([myfield],'Interpreter','None')
title('individual schnitz traces')
grid on

%define colors
myColormap=[1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0]; 
number_of_colors=size(myColormap,1);

%loop over all schnitzes
for cell=1:length(myschnitzcells)
    curX=myschnitzcells(cell).(myframes);
    curY=myschnitzcells(cell).(myfield);
    curNumRep=zeros(size(curX))+cell;  % vector in which the current schnitznumber stands. length = #frames it exists
    if ~isnan(curX) & ~isnan(curY) & length(curX)==length(curY)
        %determine color
        coloridx=mod((cell-1),(number_of_colors-1))+1;
        %plot
        plot (curX, curY,'.-','Color',myColormap(coloridx,:))
        %collect data
        alldata=[alldata; curNumRep', curX', curY'];
    end
end

%get average time trace for all schnitzes
% no fancy weighing!
framesunique=unique(alldata(:,2));
myfieldmeantrace=zeros(size(framesunique));
for step=1:length(framesunique)
    % find data points for this frame number
    idx=find(alldata(:,2)==framesunique(step));
    submyfield=alldata(idx,3);
    myfieldmeantrace(step)=mean(submyfield);
end

plot(framesunique,myfieldmeantrace,'.-k','Linewidth',2,'MarkerSize',15)


%save image
figureFileName=['TimeTraces_' myframes '_vs_' myfield];
saveSameSize(fig1,'file',[p.DJK_saveDir figureFileName '.png'], 'format', 'png');


