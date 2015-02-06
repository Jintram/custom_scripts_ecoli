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
chosenStrain = 's732';
chosenRep = 'r2';
current_branchData = myPhosphoData.(chosenStrain).(chosenRep).branchData;
muFieldName = myPhosphoData.(chosenStrain).(chosenRep).theMuField;
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
plot(current_branchData(branch_nr).frame_nrs, current_branchData(branch_nr).(muFieldName), '-o', 'color',preferredcolors(2,:),'LineWidth',1);
plot(current_branchData(branch_nr).frame_nrs, current_branchData(branch_nr).length_fitNew, '-o', 'color',preferredcolors(3,:),'LineWidth',1);


%% length normal or semilog 
h = figure(4);
set(h, 'Position', [100 100 800+100 600+100]); clf; 
grid on;

% main data
fr = current_branchData(branch_nr).frame_nrs;
ti = current_branchData(branch_nr).time;
lengths = current_branchData(branch_nr).length_fitNew;
rates = current_branchData(branch_nr).(muFieldName);

%[ax,hline1,hline2] = plotyy(fr,lengths,fr,rates,'plot','plot');
[ax,hline1,hline2] = plotyy(ti,lengths,ti,rates,'plot','plot');

hold on; % note this has to be AFTER first plots!

axes(ax(2));
lh = line([min(ti),max(ti)],[0,0],'LineWidth',2,'Color',[0 0 0])

set(gca,'FontSize',20);                           
set(ax,'FontSize',20)

set(hline1,'LineStyle','x','Color',preferredcolors(2,:),'LineWidth',2,'MarkerSize',10);
set(hline2,'LineStyle','+','Color',preferredcolors(3,:),'LineWidth',2,'MarkerSize',10);
set(ax(1),'YColor',preferredcolors(2,:))
set(ax(2),'YColor',preferredcolors(3,:))
                           
if 0, title('Length and growth rate'); end;
ylabel(ax(1),'Length ({\mu}m)','Color',preferredcolors(2,:));
ylabel(ax(2),'Growth rate (doublings/hr)','FontSize',20,'Color',preferredcolors(3,:));

ylim(ax(1),[-.5,5]);
ylim(ax(2),[-.5,5]);
xlim(ax(1),[min(ti),max(ti)]);
xlim(ax(2),[min(ti),max(ti)]);

xlhand = get(gca,'xlabel')
set(xlhand,'string','Time (mins)','fontsize',20)

legend([hline1 hline2], 'Length','Growth rate')

% Plot schnitz info also
if 0
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
else
    disp('Remember there''s the option of plotting schnitznumbers!');
end

%% same as above but suplots
FONTSIZE = 20;

h = figure(4);
set(h, 'Position', [100 100 800+100 600+100]); clf; 
grid on;

% main data
fr = current_branchData(branch_nr).frame_nrs;
ti = current_branchData(branch_nr).time;
lengths = current_branchData(branch_nr).length_fitNew;
rates = current_branchData(branch_nr).(muFieldName);



hold on; % note this has to be AFTER first plots!

set(gca,'FontSize',20);                           

% Plot lengths
subplot(2,1,1);
hline1 = plot(ti,lengths,'.');
set(hline1,'Color',preferredcolors(1,:),'LineWidth',2,'MarkerSize',10);
set(gca, 'FontSize', FONTSIZE);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% set(gca, 'xgrid', 'on', 'GridLineStyle', '-'); % ,'Xcolor',[.7 .7 .7]
xlim([min(ti),max(ti)]);
ylim([min(lengths)*.95,max(lengths)*1.05]);
ylabel('Length ({\mu}m)','FontSize',20);

% Plot growth rates
subplot(2,1,2);
hline2 = plot(ti,rates,'.-');
set(hline2,'Color',preferredcolors(1,:),'LineWidth',2,'MarkerSize',10);
set(gca, 'FontSize', FONTSIZE);
% set(gca, 'xgrid', 'on', 'GridLineStyle', '-');
xlim([min(ti),max(ti)]);
ylim([0,max(rates)*1.05]);
ylabel('Growth rate (doublings/hr)','FontSize',FONTSIZE);
xlabel('Time (mins)','fontsize',FONTSIZE);

%xlhand = get(gca,'xlabel')
%set(xlhand,'string','Time (mins)')

% Plot schnitz info also
if 1
    % obtain framenrs were schnitznrs suddenly change
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
        text(time(schnitzswitchframes(idx)), lengths(schnitzswitchframes(idx)), [num2str(uniq_schn(idx))],'FontSize',FONTSIZE);%,'Parent',ax(1));

        % % Frame nrs
        % text(schnitzswitchframes(idx), mypos-eps,...
        % ['fr=',num2str(schnitzswitchframes(idx))],'FontSize',15); 
        
        % For later use
        divisiontimes(end+1) = time(schnitzswitchframes(idx));
    end
else
    disp('Remember there''s the option of plotting schnitznumbers!');
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
    muFieldName = myPhosphoData.(chosenStrain).(chosenRep).theMuField;
    plot(current_branchData(branch_nr).time, current_branchData(branch_nr).(muFieldName), '-', 'color',[0.5,0.5,0.5],'LineWidth',2);
    %plot(current_branchData(branch_nr).time, current_branchData(branch_nr).length_fitNew, '-o', 'color',[0.5,0.5,0.5],'LineWidth',1);
end

plot([0,max(current_branchData(highlight_branch_nr).time)],[0,0],'k-')
plot(current_branchData(highlight_branch_nr).time, current_branchData(highlight_branch_nr).(muFieldName), '-', 'color',preferredcolors(2,:),'LineWidth',3);
%plot(current_branchData(highlight_branch_nr).time, current_branchData(highlight_branch_nr).length_fitNew, '-o', 'color',preferredcolors(2,:),'LineWidth',1);

% Plot division times
plot(divisiontimes,0,'^','MarkerFaceColor',preferredcolors(2,:),'MarkerEdgeColor','k','MarkerSize',12);

ylim([0,2]);
xlim([0,max(current_branchData(highlight_branch_nr).time)]);

print(h,'-dtiff','-r600','C:\Users\wehrens\Desktop\figure.tif')

% semilogy([0,max(branchData(bac_nr).frame_nrs)],[0,0],'k-')
% semilogy(branchData(bac_nr).frame_nrs, branchData(bac_nr).muP11_all, '-o', 'color',preferredcolors(1,:),'LineWidth',1);
% semilogy(branchData(bac_nr).frame_nrs, branchData(bac_nr).length_fitNew, '-o', 'color',preferredcolors(2,:),'LineWidth',1);




%% Plot multiple branches, highlight one, dual plot now
%===
current_branchData = myPhosphoData.('s732').('r1').branchData;
highlight_branch_nr = 150;
% branch 50 shows weird behavior!


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

% Plotting all branches in gray
h = figure(4); clf;
set(h, 'Position', [100 100 800+100 600+100]); 
grid on;

for branch_nr = 1:length(myPhosphoData.('s732').('r1').branchData)
    muFieldName = myPhosphoData.(chosenStrain).(chosenRep).theMuField;
    
    % obtain data for branch
    fr = current_branchData(branch_nr).frame_nrs;
    ti = current_branchData(branch_nr).time;
    lengths = current_branchData(branch_nr).length_fitNew;
    rates = current_branchData(branch_nr).(muFieldName);
    
    set(gca,'FontSize',20);                           

    % Plot lengths
    subplot(2,1,1);
    hline1 = plot(ti,lengths,'.');
    set(hline1,'Color',[.8 .8 .8],'LineWidth',1,'MarkerSize',10);
    hold on; % note this has to be AFTER first plots!
    
    % Plot growth rates
    subplot(2,1,2);
    hline2 = plot(ti,rates,'.-');
    set(hline2,'Color',[.8 .8 .8],'LineWidth',1,'MarkerSize',10);
    hold on; % note this has to be AFTER first plots!
    
    %plot(current_branchData(branch_nr).time, current_branchData(branch_nr).length_fitNew, '-o', 'color',[0.5,0.5,0.5],'LineWidth',1);
end

% obtain data for higlighted branch
fr = current_branchData(highlight_branch_nr).frame_nrs;
ti = current_branchData(highlight_branch_nr).time;
lengths = current_branchData(highlight_branch_nr).length_fitNew;
rates = current_branchData(highlight_branch_nr).(muFieldName);

% Plot highlighted length line, 
subplot(2,1,1);
hline1 = plot(ti,lengths,'.');
set(hline1,'Color',preferredcolors(1,:),'LineWidth',2,'MarkerSize',10);
hold on; % note this has to be AFTER first plots!
% Set some stuff for highlighted line
set(gca, 'FontSize', FONTSIZE);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% set(gca, 'xgrid', 'on', 'GridLineStyle', '-'); % ,'Xcolor',[.7 .7 .7]
xlim([min(ti),max(ti)]);
ylim([min(lengths)*.95,max(lengths)*1.05]);
ylabel('Length ({\mu}m)','FontSize',20);

% Plot highlighted growth rate line, set some stuff for this plot
subplot(2,1,2);
hline2 = plot(ti,rates,'-');
set(hline2,'Color',preferredcolors(1,:),'LineWidth',2,'MarkerSize',10);
hold on; % note this has to be AFTER first plots!
% Set some stuff for highlighted line
set(gca, 'FontSize', FONTSIZE);
% set(gca, 'xgrid', 'on', 'GridLineStyle', '-');
xlim([min(ti),max(ti)]);
ylim([0,max(rates)*1.05]);
ylabel('Growth rate (doublings/hr)','FontSize',FONTSIZE);
xlabel('Time (mins)','fontsize',FONTSIZE);


%plot([0,max(current_branchData(highlight_branch_nr).time)],[0,0],'k-')
%plot(current_branchData(highlight_branch_nr).time, current_branchData(highlight_branch_nr).(muFieldName), '-', 'color',preferredcolors(2,:),'LineWidth',3);
%plot(current_branchData(highlight_branch_nr).time, current_branchData(highlight_branch_nr).length_fitNew, '-o', 'color',preferredcolors(2,:),'LineWidth',1);

% Plot division times
plot(divisiontimes,0,'^','MarkerFaceColor',preferredcolors(1,:),'MarkerEdgeColor','k','MarkerSize',12);

ylim([0,2]);
xlim([0,max(current_branchData(highlight_branch_nr).time)]);

print(h,'-dtiff','-r600','C:\Users\wehrens\Desktop\figure.tif')

% semilogy([0,max(branchData(bac_nr).frame_nrs)],[0,0],'k-')
% semilogy(branchData(bac_nr).frame_nrs, branchData(bac_nr).muP11_all, '-o', 'color',preferredcolors(1,:),'LineWidth',1);
% semilogy(branchData(bac_nr).frame_nrs, branchData(bac_nr).length_fitNew, '-o', 'color',preferredcolors(2,:),'LineWidth',1);


%% 

% For autocorrelations, see phospho_autocorrelation.m.






