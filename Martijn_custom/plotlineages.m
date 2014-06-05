% Script originally by NW to plot mu(t) for multiple branches in one plot.
% (Obtained from NW at 2014/06/05)


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
    lineage=trimmed_branches(i);
    %plot(lineage.C_time,lineage.noise_dY5_sum_dt,'-k','LineWidth',1)
    plot(lineage.dY5_time,lineage.dY5,'-k','LineWidth',1)
    xlabel('time [min]','FontSize',12)
    ylabel('noise YFP production rate [a.u.]','FontSize',12)
    set(gca,'xlim',timelim)
    %ylim([0 6000])
end

% plot highlighted lineage
if ~isempty(specline)
    %lineage=trimmed_branches(specline);
    lineageRate=trimmed_branches(specline);
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