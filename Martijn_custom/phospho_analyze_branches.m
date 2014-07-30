% phoshpo_analyze_branches


%% Loop over all loaded data and make branch structures.
% ===
disp('Starting..');

% loop over my datastruct
% first per strain
strains = fieldnames(myPhosphoData);
for str_idx = 1:length(strains)

    currentstrain = char(strains(str_idx));   
        
    % then per replicate
    reps = fieldnames(myPhosphoData.(currentstrain));
    for rep_idx = 1:length(reps)
        
        currentrep = char(reps(rep_idx));
        
        % calculate branch data for each strain, rep
        myPhosphoData = phospho_getbranchdata(myPhosphoData,currentstrain,currentrep);
        
        disp([currentstrain, ', rep ',currentrep,' done ..'])
        
    end   
    
end

disp('All branchData done!')


%% choose dataset to plot
%===
current_branchData = myPhosphoData.('s732').('r1').branchData;
branch_nr = 50;

% Data where there are negative values:
%{
current_branchData = myPhosphoData.('s732').('r2').branchData;
bac_nr = 150;
%}

% plot 
%===

%{
figure(1); clf; hold on;
plot(current_branchData(bac_nr).time, current_branchData(bac_nr).muP11_all, '-o', 'color',preferredcolors(1,:),'LineWidth',1);

figure(2); clf; hold on;
plot(current_branchData(bac_nr).time, current_branchData(bac_nr).length_fitNew, '-o', 'color',preferredcolors(1,:),'LineWidth',1);
%}

% normal scale 
figure(3); clf; hold on;
set(gca,'FontSize',20);
title('Length and speed');
xlabel('Frames');
ylabel('Length (um)');
ylabel('Speed (doublings/hr)');
plot([0,max(current_branchData(branch_nr).frame_nrs)],[0,0],'k-')
plot(current_branchData(branch_nr).frame_nrs, current_branchData(branch_nr).muP11_all, '-o', 'color',preferredcolors(1,:),'LineWidth',1);
plot(current_branchData(branch_nr).frame_nrs, current_branchData(branch_nr).length_fitNew, '-o', 'color',preferredcolors(2,:),'LineWidth',1);



%% length semilog 
figure(4); clf; 

% main data
fr = current_branchData(branch_nr).frame_nrs;
ti = current_branchData(branch_nr).time;
lengths = current_branchData(branch_nr).length_fitNew;
rates = current_branchData(branch_nr).muP11_all;

%[ax,hline1,hline2] = plotyy(fr,lengths,fr,rates,'plot','plot');
[ax,hline1,hline2] = plotyy(ti,lengths,ti,rates,'plot','plot');

hold on; % note this has to be AFTER first plots!

axes(ax(2));
lh = line([min(ti),max(ti)],[0,0])
set(lh,'Color',preferredcolors(2,:))

set(gca,'FontSize',20);                           
set(ax,'FontSize',20)

set(hline1,'LineStyle','-','Marker','o','Color',preferredcolors(1,:),'LineWidth',2);
set(hline2,'LineStyle','-','Marker','o','Color',preferredcolors(2,:),'LineWidth',2);
set(ax(1),'YColor',preferredcolors(1,:))
set(ax(2),'YColor',preferredcolors(2,:))
                           
if 0, title('Length and speed'); end;
ylabel(ax(1),'Length ({\mu}m)','Color',preferredcolors(1,:));
ylabel(ax(2),'Speed (doublings/hr)','FontSize',20,'Color',preferredcolors(2,:));

ylim(ax(1),[-.5,5]);
ylim(ax(2),[-.5,5]);
xlim(ax(1),[min(ti),max(ti)]);
xlim(ax(2),[min(ti),max(ti)]);

xlhand = get(gca,'xlabel')
set(xlhand,'string','Time (mins)','fontsize',20)

% Plot schnitz info also
if 0
    [uniq_schn,ia]=unique(current_branchData(branch_nr).schnitzNrs);
    schnitzswitchframes = fr(ia);
    
    mypos=4; 
    % % add line to give schnitzes
    % lh = line(schnitzswitchframes,ones(1,length(uniq_schn))*mypos);
    % set(lh, 'Marker','o');
    dx = 5;
    for idx = [1:length(uniq_schn)]
        % % fixed location
        % text(schnitzswitchframes(idx), mypos+eps, [num2str(uniq_schn(idx))],'FontSize',15);
        % at y-value
        text(schnitzswitchframes(idx)+dx, lengths(schnitzswitchframes(idx)), [num2str(uniq_schn(idx))],'FontSize',15,'Parent',ax(1));

        % % Frame nrs
        % text(schnitzswitchframes(idx), mypos-eps,...
        % ['fr=',num2str(schnitzswitchframes(idx))],'FontSize',15); 
    end
end


% semilogy([0,max(branchData(bac_nr).frame_nrs)],[0,0],'k-')
% semilogy(branchData(bac_nr).frame_nrs, branchData(bac_nr).muP11_all, '-o', 'color',preferredcolors(1,:),'LineWidth',1);
% semilogy(branchData(bac_nr).frame_nrs, branchData(bac_nr).length_fitNew, '-o', 'color',preferredcolors(2,:),'LineWidth',1);

%% Continue to plot crosscorrelations

branch_groups_test = {myPhosphoData.('s735').('r1').branchData,myPhosphoData.('s735').('r2').branchData};
DJK_plot_crosscorrelation_standard_error_store(myPhosphoData.('s735').('r1').p,branch_groups_test,'muP11_all', 'muP11_all','selectionName','mytestname','timeField','time','onScreen',1);








% After we have data available in separate branches, we want to
% calculate auto correlation.
%{
branches = DJK_addToBranches_noise(p, branchData,'dataFields',{'dR5_time'  'R_time'  'muP11_fitNew_atdR5' 'muP11_fitNew_atdR5_cycCor' 'dR5_cycCor'  'dG5_cycCor' 'dG5' 'dR5' });
trimmed_branches = DJK_trim_branch_data(branches,4);
branch_groups = DJK_divide_branch_data(trimmed_branches);

DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dG5_cycCor', 'noise_muP11_fitNew_atdR5_cycCor','selectionName',name_rm_branch,'timeField','R_time');
%}






