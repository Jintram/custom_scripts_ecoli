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
current_branchData = myPhosphoData.('s732').('r2').branchData;
branch_nr = 150;

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
if 1
    [uniq_schn,ia]=unique(current_branchData(branch_nr).schnitzNrs);
    time=current_branchData(branch_nr).time;
    schnitzswitchframes = fr(ia);
    
    mypos=4; 
    % % add line to give schnitzes
    % lh = line(schnitzswitchframes,ones(1,length(uniq_schn))*mypos);
    % set(lh, 'Marker','o');
    dx = 5;
    divisiontimes = []; % to store times of division for later use
    for idx = [1:length(uniq_schn)]
        % % fixed location
        % text(schnitzswitchframes(idx), mypos+eps, [num2str(uniq_schn(idx))],'FontSize',15);
        % at y-value
        text(time(schnitzswitchframes(idx)), lengths(schnitzswitchframes(idx)), [num2str(uniq_schn(idx))],'FontSize',15,'Parent',ax(1));

        % % Frame nrs
        % text(schnitzswitchframes(idx), mypos-eps,...
        % ['fr=',num2str(schnitzswitchframes(idx))],'FontSize',15); 
        
        % For later use
        divisiontimes(end+1) = time(schnitzswitchframes(idx));
    end
end



%% Plot multiple branches, highlight one
%===
current_branchData = myPhosphoData.('s732').('r1').branchData;
highlight_branch_nr = 50;

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
h=figure(3); clf; hold on;
set(gca,'FontSize',20);
%title('WT Colony growth curves');
xlabel('Time');
ylabel('Speed (doublings/hr)');

for branch_nr = 1:length(myPhosphoData.('s732').('r1').branchData)
    plot(current_branchData(branch_nr).time, current_branchData(branch_nr).muP11_all, '-', 'color',[0.5,0.5,0.5],'LineWidth',2);
    %plot(current_branchData(branch_nr).time, current_branchData(branch_nr).length_fitNew, '-o', 'color',[0.5,0.5,0.5],'LineWidth',1);
end

plot([0,max(current_branchData(highlight_branch_nr).time)],[0,0],'k-')
plot(current_branchData(highlight_branch_nr).time, current_branchData(highlight_branch_nr).muP11_all, '-', 'color',preferredcolors(2,:),'LineWidth',3);
%plot(current_branchData(highlight_branch_nr).time, current_branchData(highlight_branch_nr).length_fitNew, '-o', 'color',preferredcolors(2,:),'LineWidth',1);

% Plot division times
plot(divisiontimes,0,'^','MarkerFaceColor',preferredcolors(2,:),'MarkerEdgeColor','k','MarkerSize',12);

ylim([0,2]);
xlim([0,max(current_branchData(highlight_branch_nr).time)]);

print(h,'-dtiff','-r600','C:\Users\wehrens\Desktop\figure.tif')

% semilogy([0,max(branchData(bac_nr).frame_nrs)],[0,0],'k-')
% semilogy(branchData(bac_nr).frame_nrs, branchData(bac_nr).muP11_all, '-o', 'color',preferredcolors(1,:),'LineWidth',1);
% semilogy(branchData(bac_nr).frame_nrs, branchData(bac_nr).length_fitNew, '-o', 'color',preferredcolors(2,:),'LineWidth',1);

%% 

% For autocorrelations, see phospho_autocorrelation.m.






