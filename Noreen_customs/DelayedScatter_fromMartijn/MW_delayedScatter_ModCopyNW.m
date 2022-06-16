% ************************************************************************
% This script plots a delayed scatter plot for every time delay in the
% cross-correlation
%
% Reason: Check if the scatter is simply linear (e.g. gaussian noise with
% correlation) "the higher the better". Or, if an optimum concentration etc.
% exists, e.g. "growth peaks for certain conc"
%    /               __
%   /        vs     /  \           (e.g. growth rate (y-axis), conc (x-axis))
%  /               /    \        
%
% The pro. most interesting delay is the peak of the cross-corr. The simple
% correlation coefficient cannot distinguish between the different shapes
%
% Note: weighing is partially implemented: a data point is allows to be
% used in REDUNCANCYALLOWED (typically=4) lineages, then it is ignored from then on
%
% ************************************************************************
% Script is written by Martijn Wehrens (2015-04)
% This is the working copy of Noreen Walker (slightly modified 2015-09)
% ************************************************************************



%% ***************************************************************************
% (1)Define paths & variables. Get data [ADJUST MANUALLY]
% ****************************************************************************

% ------------------ Manual input begin -----------------
% Run schnitzcells analysis .xlsx file to get:
% p, s_rm, fitTime      (fitTime can be changed below +-2)
%
myID = 'PrrnGFP_forThesis';  % strain description
%associatedFieldNames =  {'C_time','C6_mean_cycCor', 'muP9_fitNew_cycCor'}; 
%associatedFieldNames =  {'G_time','G6_mean_cycCor', 'muP15_fitNew_cycCor'}; 
associatedFieldNames =  {'dG5_time' 'dG5_cycCor' 'muP15_fitNew_atdR5_cycCor'};
                                % fields for crosscorr & scatter. Timefield first!!
                                % e.g. {'G_time','G6_mean_cycCor', 'muP15_fitNew_cycCor'}
YFIELDBRANCHPLOT = 2; % field plotted on y-axis of time traces (use 2 or 3)
% ------------------ Manual input end -----------------

% Extra parameters from MW script. Adjust if really wanted:
p.NW_saveDir=[ p.analysisDir 'DelayedScatterGFP_forThesis' filesep];
myOutputFolder=p.NW_saveDir;  % blubb some old naming relics
myTitle = myID; % Figure title (could be a different string)
badSchnitzes = []; % remove these schnitzes (should already be set in s_rm)
addSlowOnes = 0;            % Automatically add slow schnitzes to bad schnitzes
alreadyRemovedInMatFile=1;  % obsolete. If s_rm is your input, it's not necessary to generate it, then set this to 1
%blubb makeDtime = 0;              % !! If 1, re-calculates schnitzcells.dX_time field !!
                                    % if set to 1, extra MWDJK_dX_time filed is recalculated.
PLOTSCATTER=1;              % If zero, only cross corrs are calculated, not delayed scatter plots
PERFORMSOMECHECKS = 0; % NW: no idea what that is
PLOTLINEAGESHIST= 0;    % time traces of lineages & histogram. see cell (2)



% Make sure that NW_saveDir directory exists
if exist(p.NW_saveDir)~=7
  [status,msg,id] = mymkdir([p.NW_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.NW_saveDir ' : ' msg]);
    return;
  end
end
% Check for input errors
if ~exist('associatedFieldNames') | ~exist('myTitle') | ~exist('p') | ~exist('badSchnitzes')
    error('input not supplied.')
end
    
name_rm = 'rm'; name_all = 'all'; % necessary?


%% ***************************************************************************
% (2)Get branches/lineages. Plot time traces and histogram
% ****************************************************************************

% HIGHLIGHTSUSPICOUS = 0; not working currently
% fitTime = fitTime + [2 -2]; blubb. if changed here, make sure to reload from xlsx file

% obtain correct -1 shifted frame nrs
s_rm = MW_calculateframe_nrs_ModCopyNW(s_rm); % shifted frames (used in MW schnitzcells 
            % version, not in NW schnitzcells version yet. (backwards compatibility fix )

% obtain branches
branchData = DJK_getBranches(p,s_rm,'dataFields',{associatedFieldNames{1}, associatedFieldNames{2}, ...
    associatedFieldNames{3} }, 'fitTime', fitTime); 
name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2)) '_CrossCorrs'];

% colors for plotting branches
distinguishableColors = distinguishable_colors_CopyNW(numel(branchData)+1,[1 1 1]); 


% --------------------------------
% Plot all branches (time traces of lineages). field(YFIELDTOPLOT) vs field1 (e.g. conc vs time)
% --------------------------------
if PLOTLINEAGESHIST
    
    figure(1); clf; hold on;
    numelBranches = numel(branchData);
    for branchIdx = 1:numelBranches
        l = plot(branchData(branchIdx).(associatedFieldNames{1}), branchData(branchIdx).(associatedFieldNames{YFIELDBRANCHPLOT}),'-o','Color',distinguishableColors(branchIdx,:));
        set(l, 'LineWidth', (numelBranches-branchIdx+1)/numelBranches*10);
    end

    % make figure pretty
    xlabel(associatedFieldNames{1},'Interpreter', 'None'), ylabel(associatedFieldNames{YFIELDBRANCHPLOT},'Interpreter', 'None')
    myXlimFig1 = max(branchData(branchIdx).(associatedFieldNames{1}));
    xlim([0, myXlimFig1]);
    myYlimFig1 = [min([branchData.(associatedFieldNames{YFIELDBRANCHPLOT})]),...
                  max([branchData.(associatedFieldNames{YFIELDBRANCHPLOT})])];
    ylim([0, myYlimFig1(2)*1.5]);
    MW_makeplotlookbetter_CopyNW(12); % set font size

    % --------------------------------
    % Plot Histogram of field YFIELDTOPLOT (e.g. conc)
    % --------------------------------
    figure(2), clf, hold on
    allYdata = [branchData.(associatedFieldNames{YFIELDBRANCHPLOT})];
    [nelements, centers] = hist(allYdata,20);
    deltaY = centers(2)-centers(1);
    totalCount = numel(allYdata);
    %nelements=nelements./deltaY;
    plot(centers,nelements,'or','LineWidth',2)
    % Fit distribution
    pd=fitdist(allYdata', 'Normal');
    fittedDistrX = [min(allYdata):(max(allYdata)-min(allYdata))/100:max(allYdata)];
    fittedDistrYnorm = normpdf(fittedDistrX,pd.mu,pd.sigma);
    fittedDistrY = fittedDistrYnorm.*totalCount.*deltaY;
    plot(fittedDistrX,fittedDistrY,'k', 'LineWidth', 3);
    %probplot(allYdata)
    MW_makeplotlookbetter_CopyNW(12);

    title('CAREFUL: WIEGHING INCORRECT! FIRST CELL WEIGHED ~500-times') % NW2015-09
    xlabel('Fluor strength s (a.u.)');
    ylabel('P(s) * N * \Deltas (counts)');

    % NW: NOT WORKING..start...
    %{ 
    %  TODO make 99% confidence and plot in previous figure.
     %confidence = paramci(pd,'Alpha',.01);
     sigma2 = pd.mu + 4.*[-pd.sigma, pd.sigma];
    plot([sigma2(1),sigma2(1)],[0,myYlim],'--','Color',[.5 .5 .5],'LineWidth', 2)
    plot([sigma2(2),sigma2(2)],[0,myYlim],'--','Color',[.5 .5 .5],'LineWidth', 2)
    sigma4 = pd.mu + 4.*[-pd.sigma, pd.sigma];
    plot([sigma5(1),sigma5(1)],[0,myYlim],':','Color',[.5 .5 .5],'LineWidth', 2)
    plot([sigma5(2),sigma5(2)],[0,myYlim],':','Color',[.5 .5 .5],'LineWidth', 2)

    myYlim = max(nelements)*1.1;
    ylim([0,myYlim]);
    % NW: blubb ??? xlim([myYlimFig1(1), myYlimFig1(2)*1.5]);

    % Now also plot confidence intervals in previous figure
    figure(1), hold on;
    plot([0,myXlimFig1],[sigma2(1),sigma2(1)],'--','Color',[.5 .5 .5],'LineWidth', 2)
    plot([0,myXlimFig1],[sigma2(2),sigma2(2)],'--','Color',[.5 .5 .5],'LineWidth', 2)
    plot([0,myXlimFig1],[sigma5(1),sigma5(1)],':','Color',[.5 .5 .5],'LineWidth', 2)
    plot([0,myXlimFig1],[sigma5(2),sigma5(2)],':','Color',[.5 .5 .5],'LineWidth', 2)

    % Now list schnitzes that have suspiciously high signal:
    suspiciousBranches = []; suspiciousSchnitzes = [];
    for branchIdx = 1:numel(branchData)
        if any(branchData(branchIdx).(associatedFieldNames{YFIELDBRANCHPLOT}) > sigma5(2))
            suspiciousBranches(end+1) = branchIdx;
            if HIGHLIGHTSUSPICOUS % plot if desired
                plot(branchData(branchIdx).(associatedFieldNames{1}), branchData(branchIdx).(associatedFieldNames{YFIELDBRANCHPLOT}),'-or','LineWidth',3)
            end
            %plot(branchData(branchIdx).(associatedFieldNames{1}),
            %branchData(branchIdx).(associatedFieldNames{YFIELDBRANCHPLOT}),'-o','LineWidth',3,'Color',mycolors(c)) % MW debug

            % Find out which schnitzes
            locationsInThisBranch = find(branchData(branchIdx).(associatedFieldNames{YFIELDBRANCHPLOT}) > sigma5(2));
            suspiciousSchnitzes = [suspiciousSchnitzes branchData(branchIdx).schnitzNrs(locationsInThisBranch)];
        end
    end
    suspiciousSchnitzes = unique(suspiciousSchnitzes)
    suspiciousBranches
    %} 
    % NW: NOT WORKING ... end...
end

%% ***************************************************************************
% (3)MAIN PART. Get cross-correlations. Get delayed scatter plots [Manual
% Adjustment Optional]
% ****************************************************************************

% ----------- Adjust if wanted start ----------------------
NRCONTOURLINES = 3;       % # countour lines of equal probability in scatter plot (kernel density estimate)
SHOWPLUSMINFROMZERO = 6; % max time delay. unit: # fluo frames
% how many times can one datapoint be reused for different lineages
REDUNDANCYALLOWED = 2^2; 
% ----------- Adjust if wanted end ----------------------

% -------------------------------------------------------------
% get standard cross-correlations (and store in delayedscatter folder)
% -------------------------------------------------------------
% get branchgroups
branches = DJK_addToBranches_noise(p, branchData,'dataFields',{associatedFieldNames{1},associatedFieldNames{2},associatedFieldNames{3}});
trimmed_branches = DJK_trim_branch_data(branches,4); % divide into 4 branchgroups NW2015-09
branch_groups = DJK_divide_branch_data(trimmed_branches);

% Colony average mean has already been substracted so theoretically extra
% normalization shouldn't have an effect.
p.extraNorm=0; % NW: ???

% THIS MIGHT FAIL BECAUSE TIME IS NOT SET CORRECTLY (SHOULD BE time_at_..)
% To calculate cross correlations, additional normalization is usually
% performed, namely to filter out colony average behavior.
[CorrData,composite_corr] = DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, ['noise_' associatedFieldNames{1,2}],['noise_' associatedFieldNames{1,3}] ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',1,'DJK_saveDir',p.NW_saveDir); 



% -------------------------------------------------------------
% get delayed scatter data
% -------------------------------------------------------------
% Do we want to filter out colony average behavior for the "delayed
% scatter" plots also? Maybe do this with noise fields?
% But let's try with "raw" data first..
p.timeField = associatedFieldNames{1,1};

[dataPairsPerTau, iTausCalculated, originColorPerTau, correlationsPerTau] = ...
    MW_getdelayedscatter_ModCopyNW(p, branches, ['noise_' associatedFieldNames{1,2}], ['noise_' associatedFieldNames{1,3}], REDUNDANCYALLOWED);

% -------------------------------------------------------------
% some Extra Stuff CCs comparisons for different calculation methods (not used NW)
% -------------------------------------------------------------
%{
Plot "raw" cross cor I calculate (MW)

myfig=figure(99),clf,hold on;
l=plot(iTausCalculated,correlationsPerTau,'o-r','LineWidth',2)

% Compare two cross-corrs (DJK & MW)

myfig=figure(5),clf,hold on;
errorbar(CorrData(:,1),CorrData(:,2),CorrData(:,3),'x-','Color', [.5,.5,.5], 'LineWidth',2)
l1=plot(CorrData(:,1),CorrData(:,2),'x-k','LineWidth',2)

l2=plot(CorrData(:,1),correlationsPerTau,'o-r','LineWidth',2)

% If you recalculate correlations again w. different params, this allows
% plotting of extra line.
%l3=plot(CorrData(:,1),correlationsPerTau100,'o-','Color',[1 .5 0],'LineWidth',2)

myxlimvalues=[min(CorrData(:,1)), max(CorrData(:,1))];
xlim(myxlimvalues);
ylim([-1,1]);
plot(myxlimvalues,[0,0],'k-');

legend([l1,l2],{'DJK','MW'})

title(['DJK vs. MW R -- ' myID '_' p. movieDate  '_' p.movieName], 'Interpreter', 'none');
xlabel('\tau (hrs)');
ylabel(['R(' associatedFieldNames{1,2} ', growth) (normalized)'], 'Interpreter', 'none');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

plot([0,0],[-1,1],'-k');

saveas(myfig,[myOutputFolder 'crosscorrs_' associatedFieldNames{1,2} '.tiff']);
%}

% -------------------------------------------------------------
% Figure delayed scatters
% -------------------------------------------------------------
if ~exist('PLOTSCATTER','var'), PLOTSCATTER=1; end;

if PLOTSCATTER
    middlePosition = [ceil(numel(iTausCalculated)/2)-SHOWPLUSMINFROMZERO:ceil(numel(iTausCalculated)/2)+SHOWPLUSMINFROMZERO];
    rangeiTausCalculated = middlePosition;
    % delayIdx = 11; % 39 is middle

    myColorMap = colormap(winter(numel(iTausCalculated)));

    
    hFig = figure(3); clf; hold on;

    densities=[]; Xs=[]; Ys=[]; Zs=[];
    for delayIdx = rangeiTausCalculated

        % Rename data more convenient
        data = [dataPairsPerTau{delayIdx}(:,1), dataPairsPerTau{delayIdx}(:,2)];    

        
        % plot
        set(0,'CurrentFigure',3); clf, hold on;
        offset=100; width1=500; height1=500;
        set(hFig, 'Position', [(offset+width1) offset width1 height1]);     

        % scatter
        %for pointIdx = 1:numel(data(:,1))
        %    %plot(data(pointIdx,1),data(pointIdx,2),'.','Color',originColorPerTau{delayIdx}(pointIdx,:));%,'Color',distinguishableColors(myGrouping(i)+1,:),'MarkerSize',3);
        %    plot(data(pointIdx,1),data(pointIdx,2),'.','Color','r');%
        %end
            plot(data(:,1),data(:,2),'o','Color','r','MarkerSize',2);% alternative: quick and unicolor

        % contour (from kde)
        [bandwidth,density,X,Y] = kde2d(data);    
        [C, l1] = contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);

        % mean
        lineH = plot(mean(data(:,1)),mean(data(:,2)),'o','MarkerFaceColor','k','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);

        title(['D# = ' num2str(iTausCalculated(delayIdx)) ', R = ' sprintf('%0.3f',correlationsPerTau(delayIdx)), 10 ,myID '_' p. movieDate  '_' p.movieName], 'Interpreter', 'none')    

        xlabel(['Delta ' associatedFieldNames{1,2}] , 'Interpreter', 'none'); % = "noise_C6_mean_cycCor" etc
        ylabel(['Delta ' associatedFieldNames{1,3}] , 'Interpreter', 'none');
        %Set all fontsizes
        set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
        set(gca,'FontSize',15);

        midIdx = ceil(numel(iTausCalculated)/2);
        xlim([  min(dataPairsPerTau{midIdx}(:,1)), max(dataPairsPerTau{midIdx}(:,1))  ])
        ylim([  min(dataPairsPerTau{midIdx}(:,2)), max(dataPairsPerTau{midIdx}(:,2))  ])

        % create figure name (NW2015-09)
        delayinframes=str3(iTausCalculated(delayIdx));
        myfigname=['DelayScatter_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '_idx' str2(delayIdx) '_delay_' delayinframes];
        %save in different formats
        savefig(3,[p.NW_saveDir myfigname ])
        saveas(3,[p.NW_saveDir myfigname '.tiff']);
        saveas(3,[p.NW_saveDir myfigname '.pdf']);
        
    end

    figure(3);

    % [not used] average point (used for legend too)
    %{
    legendLines = []; previous = 0;
    for i = 1:numberOfDataFiles        
        lineH = plot(mean(data(:)),mean(data(2,:)),'o','MarkerFaceColor',distinguishableColors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);
    if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
    end    
    legend( legendLines, legendDescriptions,'Location','northeast');
    %}

    
    
end

% delete NW_saveDir
 rmfield(p,'NW_saveDir')

%% ************************************************************************
%  Some random checks. not used in this version (NW)
%  ************************************************************************

%{
if PERFORMSOMECHECKS
    % distribution of times 
    % equivalent of 'hist' calculated manually because times have interval,
    % which complicates use of hist function.
    fieldcount=0; R = struct;
    % loop over timefields of interest
    for myfield = {'Y_time','MWDJK_dY_time', 'time'}

        % administration
        myfield = char(myfield);
        fieldcount = fieldcount+1;

        % get all possible values 
        R(fieldcount).YValues = unique([s_rm.(myfield)]);
        % count how many of each values are in the schnitz struct
        R(fieldcount).countYValues = [];
        for value=R(fieldcount).YValues 
            R(fieldcount).countYValues(end+1) = sum([s_rm.(myfield)]==value);
        end

    end
    % plot first two fields of interest
    figure, clf; hold on;
    if numel(R)>2, plot(R(3).YValues, R(3).countYValues, 'xk'); end
    plot(R(1).YValues, R(1).countYValues, 'or','LineWidth',2);
    plot(R(2).YValues, R(2).countYValues, 'ob','LineWidth',2);
end
%}
