
% REMOVES A CELL NUMBER IN Lc THAT DOES NOT CORRESPOND TO A CELL (number
% exists, but area is empty) andrenumbers the highest cell number with this
% obsolete number
% press any key during run to proceed

%load('pos1cropseg309.mat')
%waitforbuttonpress 
% insert own matrix name (Lc) in 2nd and last line!
% *********************ADJUST*******************************
framenr=580;
obsoletecellnumber=475;       % empty cell number
% ************************************************************

name= [p.segmentationDir,p.movieName,'seg',str3(framenr)];
load(name);

mat=Lc;  % matrix name
imagesc(mat==obsoletecellnumber)
waitforbuttonpress
maxcellnumber=max(max(mat));  % highest cell number
for i=1:size(mat,1)
    for j=1:size(mat,2)
        if mat(i,j)==maxcellnumber
            mat(i,j)=obsoletecellnumber;
        end
    end
end
Lc=mat;
imagesc(mat==obsoletecellnumber);
clear i j mat maxcellnumber obsoletecellnumber;
savelist=['''phsub'',''LNsub'',''rect'',''timestamp'',''phaseFullSize'',''tempsegcorrect'',''Lc'''];
Lname = [p.movieName, 'seg', str3(framenr)];
eval(['save(''',p.segmentationDir,Lname,''',',savelist,');']);

%**************************************************************
%**************************************************************
%%
% CHECKS TIMESTAMPS OF IMAGE (seg) FILES
% control, if loop length in experiment was long enough and images are
% taken after equal time lapses

% *** adjust ***
segdirectory='U:\Experiments\2011-11-08\pos1crop\segmentation\';
%segdirectory='D:\DummyExp\2011-xx-xx\pos3crop\segmentation\';
looplength=50; % length of loop in sec
framefreq=1;   % every framefreq'th frame is segmented (every frame: 1)
% **************

imagenames=sprintf('%spos*',segdirectory);
images=dir(imagenames);
numimages=length(images);
datetime=zeros(1,numimages);
for i=1:numimages
    file=sprintf('%s%s',segdirectory,images(i).name);
    load(file,'timestamp');
    datetime(i)=timestamp;
end

% convert time to seconds: copied and adjusted from Daan
s = sort([datetime(find(datetime>0))]);
firstmoment = s(1);
time_sec=[];
for i = 1:length(datetime),
  if (datetime(i)>0)
    tdiff = datetime(i) - firstmoment;
    [y,m,d,h,mi,s] = datevec(tdiff);    
     time_sec(i) = 60*(s/60 + mi + h*60 + d*24*60);
  else
      time_sec(i)=0;
  end
end

%looplength=75;
% lost time (=0 if loop length was long enough)
deltatime=zeros(1,length(time_sec));
for i=1:length(time_sec)    
    deltatime(i)=time_sec(i)-looplength*(i-1)*framefreq;
end

% plot the time that images are behind
figure(1)
clf
plot(deltatime,'.')
xlabel('frame number');
ylabel('time too late [sec]');
hold on
grid on
hold off

% average time lost per frame
% approximativ (ueberschaetzt) wenn letzte Runde ein Fluo-Bild war
timelate_total=deltatime(length(deltatime));
timelate_perframe=timelate_total/length(deltatime);
str=sprintf('%4.4f sec behind within %2.0f frames...%4.4f sec per frame',timelate_perframe,framefreq,timelate_perframe/framefreq); disp(str);

%**************************************************************
%*******************************************************************
%%
% PLOTS SINGLE LINEAGE TIME TRACES OF YFP,CFP,mu,CFPprodRate,YFPprodRate
% ***TRIMMED_BRANCHES**  must have been acquired before (see excel file)
timelim=[0 2000];      % time range in [min]
lineagenumbers=[1:1000];%nicebranches;   % which lineages to plot
specline=[6 ];            % highlighted lineage  (use [] if no lineage should be highlighted

% plot YFP
figure(1)
clf
hold on
for run=1:length(lineagenumbers)
    i=lineagenumbers(run);
    lineage=trimmed_branches(i);
    plot(lineage.Y_time,lineage.Y6_mean_cycCor,'-k')
    xlabel('time [min]')
    ylabel('YFP [a.u.]')
    set(gca,'xlim',timelim)
end


% plot CFP
figure(2)
clf
hold on
for run=1:length(lineagenumbers)
    i=lineagenumbers(run);
    lineage=trimmed_branches(i);
    plot(lineage.Y_time,lineage.C6_mean_cycCor,'-k')
    xlabel('time [min]')
    ylabel('CFP [a.u.]')
    set(gca,'xlim',timelim)
end

%% and more lineages:
% plot area
figure(1)
clf
hold on
for run=1:length(lineagenumbers)
    i=lineagenumbers(run);
    lineage=trimmed_branches(i);
    plot(lineage.time,lineage.area,'-k')
    xlabel('time [min]')
    ylabel('area [a.u.]')
    set(gca,'xlim',timelim)
end


% plot Length
figure(2)
clf
hold on
for run=1:length(lineagenumbers)
    i=lineagenumbers(run);
    lineage=trimmed_branches(i);
    plot(lineage.time,lineage.length_fitNew,'-k')
    xlabel('time [min]')
    ylabel('length [a.u.]')
    set(gca,'xlim',timelim)
end


% plot mu
figure(3)
clf
hold on
for run=1:length(lineagenumbers)
    i=lineagenumbers(run);
    lineage=trimmed_branches(i);
    plot(lineage.C_time,lineage.muP15_fitNew,'-k','LineWidth',1)
    xlabel('time [min]','FontSize',12)
    ylabel('mu [dbl/hour]','FontSize',12)
    set(gca,'xlim',timelim)
end
%%

timelim=[250 850];      % time range in [min]
lineagenumbers=[3:5];%nicebranches;   % which lineages to plot
specline=[6 ];            % highlighted lineage  (use [] if no lineage should be highlighted


% plot YFP production
figure(4)
clf
hold on
for run=1:length(lineagenumbers)
    i=lineagenumbers(run);
    lineage=trimmed_branchesRate(i);
    %plot(lineage.C_time,lineage.noise_dY5_sum_dt,'-k','LineWidth',1)
    plot(lineage.dY5_time,lineage.dY5,'-k','LineWidth',1)
    xlabel('time [min]','FontSize',12)
    ylabel('noise YFP production rate [a.u.]','FontSize',12)
    set(gca,'xlim',timelim)
    %ylim([0 6000])
end

% % plot CFP production
%figure(5)
%clf
%hold on
%for run=1:length(lineagenumbers)
%    i=lineagenumbers(run);
%    lineage=trimmed_branchesRate(i);
%    %plot(lineage.C_time,lineage.noise_dC5_sum_dt,'-k')
%    plot(lineage.dY5_time,lineage.dC5,'-k')    
%    xlabel('time [min]')
%    ylabel('noise CFP production rate [a.u.]')
%    set(gca,'xlim',timelim)
%end


% plot highlighted lineage
if ~isempty(specline)
    %lineage=trimmed_branches(specline);
    lineageRate=trimmed_branchesRate(specline);
  %  figure(1)
  %  %plot(lineage.Y_time,lineage.Y6_mean,'-r','LineWidth',2)
  %  plot(lineage.time,lineage.area,'-r','LineWidth',2)
  %  figure(2)
  %  %plot(lineage.C_time,lineage.C6_mean,'-r','LineWidth',2)
  %  plot(lineage.time,lineage.length_fitNew,'-r','LineWidth',2)
  %  
  %  figure(3)
  %  plot(lineage.C_time,lineage.muP15_fitNew,'-r','LineWidth',2)
    figure(4)
    plot(lineage.dY5_time,lineage.dY5,'-r','LineWidth',2)
   % figure(5)
   % plot(lineage.C_time,lineage.dC5_sum_dt,'-r','LineWidth',2)

end


%%
% searches for nice looking branches that stay within certain Y-range
timelim=[300 500] ;
nicebranches=[];

for i=1:length(trimmed_branches)
    intimeframes=find(trimmed_branches(i).Y_time>timelim(1) & trimmed_branches(i).Y_time<timelim(2));
    weirdo=find(trimmed_branches(i).dY5_sum_dt(intimeframes)<0 | (trimmed_branches(i).dY5_sum_dt(intimeframes)>6000));
 
    if isempty(weirdo)
        disp([num2str(i)])
        nicebranches=[nicebranches;i];
    end
end