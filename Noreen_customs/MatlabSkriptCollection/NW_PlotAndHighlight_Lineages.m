% Miniscript to find lineages for the simultaneous-bursting at basal
% expression
% It is based on the search for the last(!) schnitz in a lineage and does
% not search for intermediate schntiznrs

% *******
% THIS IS A COPY OF HelpScriptPlotLineages.m IN THE EXTIRNSIC NOISE PROJECT
% FOLDER
% *******

% 2012-05-16pos1
% time range of the jpg images: frame431=715min, frame511=848min

% ****** SEE ALSO COMMENTS BELOW BEFORE CONCENTRATION PLOT!!! ***********
% load from \\biofysicasrv\Users2\Walker\ExtrinsicNoise\Analysis\BasalExpression_SimultaneousBursting\Data
% "BranchesCombinedConcRate.mat"
% **********************************************************************
% the .mat files in this folder contain "trimmed_branches" and
% "lastschnitzes" which is the last schnitznumber of each branch [in the
% same order as trimmed_branches!] (name indicates whether conc "C6Y6" or
% Rate "dC5dY5"

% now SELECT SCHNITZ 1113 & 1118 which shall be HIGHLIGHTED lineages [their
% parent is 535]
schn1=1113;
schn2=1118;


% CONCENTRATION
idx1conc=find(lastschnitzesC6Y6==schn1); % first special lineage
idx2conc=find(lastschnitzesC6Y6==schn2); % 2nd special lineage

% check that order has not been messed up:
if trimmed_branchesC6Y6(idx1conc).schnitzNrs(end)~=schn1 | trimmed_branchesC6Y6(idx2conc).schnitzNrs(end)~=schn2
    error('Sth wrong with schnitznr association!')
end

%RATE
idx1rate=find(lastschnitzesdC5dY5==schn1); % first special lineage
idx2rate=find(lastschnitzesdC5dY5==schn2); % 2nd special lineage

% check that order has not been messed up:
if trimmed_branchesdC5dY5(idx1rate).schnitzNrs(end)~=schn1 | trimmed_branchesdC5dY5(idx2rate).schnitzNrs(end)~=schn2
    error('Sth wrong with schnitznr association!')
end

% -----------------------------------------------------------------
% RUN FOR CONCENTRATION
% -----------------------------------------------------------------
% plot lineages & highlighted lineages

% basically copied from the 'littlehelpers' script

timelim=[500 870];      % time range in [min] % works: [500 870], images: 715-848min
lineagenumbersconc=[1:518];%[350:400];  % which lineages to plot %works: [350:400] 
speclineconc=[idx1conc idx2conc ];            % highlighted lineage [ all in red]  (use [] if no lineage should be highlighted
otherspeclineconc=[250 267];%[200 250];    % other highlighted lineages [all in different colors]

% *************** NOTE *********************
% For production rate the schnitzes corresponding to the lineages of of 'lineagnumbers' and 
% 'otherspecline' will be infered below so that the "same" lineages (=with
% the same cells) are plotted
% It automatically checks whether lineages exist for both conc&rate and
% displays warning/error
% ******************************************


if ~isempty(otherspeclineconc)
    % myColor needs to be extended if more than 8 lineages plotted
    myColor=[0 0.8 0 ; 0 0 1; 1 0.6 0.2;  0 0.8 0.8;  0.6 0 0.4; 0.8 0.3 0; 0.8 0.1 0; 0.4 0.2 0.6];
end

% plot YFP conc
figure(1)
clf
hold on
set(gcf,'WindowStyle','docked')
for run=1:length(lineagenumbersconc)
    i=lineagenumbersconc(run);
    lineage=trimmed_branchesC6Y6(i);
    plot(lineage.Y_time,lineage.Y6_mean_cycCor,'-','Color',[0.4 0.4 0.4],'LineWidth',2)
    xlabel('time [min]')
    ylabel('YFP conc [a.u.]')
    set(gca,'xlim',timelim)
        if ~isempty(otherspeclineconc)
        for ss=1:length(otherspeclineconc)
            lineage=trimmed_branchesC6Y6(otherspeclineconc(ss));
            plot(lineage.Y_time,lineage.Y6_mean_cycCor,'-','Color',myColor(ss,:),'LineWidth',3)
        end
        end
    if ~isempty(speclineconc)
        for ss=1:length(speclineconc)
            lineage=trimmed_branchesC6Y6(speclineconc(ss));
            plot(lineage.Y_time,lineage.Y6_mean_cycCor,'-r','LineWidth',3)
        end
    end
            
end


% plot CFP conc
figure(2)
clf
hold on
set(gcf,'WindowStyle','docked')
for run=1:length(lineagenumbersconc)
    i=lineagenumbersconc(run);
    lineage=trimmed_branchesC6Y6(i);
    plot(lineage.Y_time,lineage.C6_mean_cycCor,'-','Color',[0.4 0.4 0.4],'LineWidth',2)
    xlabel('time [min]')
    ylabel('CFP conc [a.u.]')
    set(gca,'xlim',timelim)
    if ~isempty(otherspeclineconc)
        for ss=1:length(otherspeclineconc)
            lineage=trimmed_branchesC6Y6(otherspeclineconc(ss));
            plot(lineage.Y_time,lineage.C6_mean_cycCor,'-','Color',myColor(ss,:),'LineWidth',3)
        end
    end
    if ~isempty(speclineconc)
        for ss=1:length(speclineconc)
            lineage=trimmed_branchesC6Y6(speclineconc(ss));
            plot(lineage.Y_time,lineage.C6_mean_cycCor,'-r','LineWidth',3)
        end
    end
end


% -----------------------------------------------------------------
% RUN FOR PRODUCTION RATES
% -----------------------------------------------------------------

% SAME AS CONC TIMELIM. timelim=[0 1000];      % time range in [min]
speclinerate=[idx1rate idx2rate ];            % highlighted lineage  (use [] if no lineage should be highlighted

% INFERED FROM CONC: lineagenumbersrate=[300:400];  % which lineages to plot
schnitzeslineagenumbers=lastschnitzesC6Y6(lineagenumbersconc);
lineagenumbersrate=find(ismember(lastschnitzesdC5dY5,schnitzeslineagenumbers));
% INFERED FROM CONC: otherspecline
schnitzesotherspecline=lastschnitzesC6Y6(otherspeclineconc);
otherspeclinerate=find(ismember(lastschnitzesdC5dY5,schnitzesotherspecline));

% check for errors
if length(lineagenumbersconc)~=length(lineagenumbersrate)
    disp('WARNING: #BACKGROUND LINEAGES ''lineagenumbers'' NOT THE SAME FOR RATE & CONC')
end
if length(otherspeclineconc)~=length(otherspeclinerate)
    error('Highlighted lineages in Conc don''t exist all for Rate!!')
end
 

% plot YFP rate
figure(5)
clf
hold on
set(gcf,'WindowStyle','docked')
for run=1:length(lineagenumbersrate)
    i=lineagenumbersrate(run);
    lineage=trimmed_branchesdC5dY5(i);
    plot(lineage.dY5_time,lineage.dY5_cycCor,'-','Color',[0.4 0.4 0.4],'LineWidth',2)
    xlabel('time [min]')
    ylabel('YFP rate [a.u.]')
    set(gca,'xlim',timelim)
    if ~isempty(otherspeclinerate)
        for ss=1:length(otherspeclinerate)
            lineage=trimmed_branchesdC5dY5(otherspeclinerate(ss));
            plot(lineage.dY5_time,lineage.dY5_cycCor,'-','Color',myColor(ss,:),'LineWidth',3)
        end
    end
    
    if ~isempty(speclinerate)
        for ss=1:length(speclinerate)
            lineage=trimmed_branchesdC5dY5(speclinerate(ss));
            plot(lineage.dY5_time,lineage.dY5_cycCor,'-r','LineWidth',3)
        end
    end        
end


% plot CFP rate
figure(6)
clf
hold on
set(gcf,'WindowStyle','docked')
for run=1:length(lineagenumbersrate)
    i=lineagenumbersrate(run);
    lineage=trimmed_branchesdC5dY5(i);
    plot(lineage.dY5_time,lineage.dC5_cycCor,'-','Color',[0.4 0.4 0.4],'LineWidth',2)
    xlabel('time [min]')
    ylabel('CFP rate [a.u.]')
    set(gca,'xlim',timelim)
    if ~isempty(otherspeclinerate)
        for ss=1:length(otherspeclinerate)
            lineage=trimmed_branchesdC5dY5(otherspeclinerate(ss));
            plot(lineage.dY5_time,lineage.dC5_cycCor,'-','Color',myColor(ss,:),'LineWidth',3)
        end
    end
    
    if ~isempty(speclinerate)
        for ss=1:length(speclinerate)
            lineage=trimmed_branchesdC5dY5(speclinerate(ss));
            plot(lineage.dY5_time,lineage.dC5_cycCor,'-r','LineWidth',3)
        end
    end
end

%%
% plot the bar for the movie pics time range
figure(1)
plot([715 848],[0 0],'-k','LineWidth',5)
figure(2)
plot([715 848],[0 0],'-k','LineWidth',5)
figure(5)
plot([715 848],[-500 -500],'-k','LineWidth',5)
figure(6)
plot([715 848],[-200 -200],'-k','LineWidth',5)
