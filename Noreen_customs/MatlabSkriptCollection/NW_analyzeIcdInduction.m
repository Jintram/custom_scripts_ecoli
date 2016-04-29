function [g,s] = NW_analyzeIcdInduction(p,schnitzcells,myname)
% plots a histogram of G6_mean (gfp icd at 100ms ilumination) and S6_mean
% (gfp icd at 1000ms illumination)
% very simple function, not foolproof!
%
% input: p, schnitzcells
%        myname: name under which histogram will be saved
%
% G6_mean HAS TO BE 100ms ILLUMINATION AND S6_mean HAS TO BE 1SEC
% ILLUMINATION!!! (otherwise titles are wrong!)


%make save directory
p.DJK_saveDir = [p.analysisDir 'schnitzcells' filesep];
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end

gfp100ms=[schnitzcells.G5_mean];
%gfp1sec=[schnitzcells.S6_mean];
%gfp1sec=[schnitzcells.D6_mean];
gfp1sec=[schnitzcells.D5_mean];
notnan=find(~isnan(gfp100ms));
gfp100ms=gfp100ms(notnan);
notnan=find(~isnan(gfp1sec));
gfp1sec=gfp1sec(notnan);
g=gfp100ms;
s=gfp1sec;

fig1=figure;
clf
hold on
hist(gfp100ms)
titlestring = sprintf('%s. 100ms illum . mean=%0.2f . std=%0.2f . std/mean=%0.2f',myname,mean(gfp100ms),std(gfp100ms),std(gfp100ms)/mean(gfp100ms));
title(titlestring,'Interpreter','None')
saveSameSize(fig1,'file',[p.DJK_saveDir myname '_0100ms.png'], 'format', 'png');
close(fig1)

fig2=figure;
clf
hold on
hist(gfp1sec)
%titlestring = sprintf('%s. 1sec illum . mean=%0.2f . std=%0.2f . std/mean=%0.2f',myname,mean(gfp1sec),std(gfp1sec),std(gfp1sec)/mean(gfp1sec));
titlestring = sprintf('%s. 60ms illum . mean=%0.2f . std=%0.2f . std/mean=%0.2f',myname,mean(gfp1sec),std(gfp1sec),std(gfp1sec)/mean(gfp1sec));
title(titlestring,'Interpreter','None')
%saveSameSize(fig2,'file',[p.DJK_saveDir myname '_1000ms.png'], 'format', 'png');
saveSameSize(fig2,'file',[p.DJK_saveDir myname '_0060ms.png'], 'format', 'png');
close(fig2)