
%% Parameters for all code ================================================

% CONFIGFILE='U:\ZZ_EXPERIMENTAL_DATA\Data_per_project\CRPcAMP\config_projectCRPcAMP.m';

PLOTOUTDIR = '\\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\A_Step5_Plots\CRPcAMP\';
EXTENSIONS = {'png','svg','fig'};

some_colors;

% =========================================================================

%% Part I =================================================================


%% Parameters for part I 
SELECTIONFIELD = 'groupID2';
UBERLISTDATASETSTOPLOT = {{'extracell800_RCRP'}, {'extracell800_s70'}}; 

UBERLISTDATASETSTOPLOT = {{'extracell800_RCRP'}, {'extracell800_s70'},{'WT_RCRP'},{'WT_s70'}}; 

TOPLOTFIELDNAMES = {{'concentrationCorrData', 'rateCorrData','concentrationDualCrossCorrData', 'rateDualCrossCorrData'},...
                    {'concentrationCorrData', 'rateCorrData'}}; 
               
ADDAVERAGETOFIRSTPLOT=1;

%%

for UberListIndex = 1:numel(UBERLISTDATASETSTOPLOT)
for DUALCOLORINDEX=2
for plotFieldIndex=1:numel(TOPLOTFIELDNAMES{DUALCOLORINDEX})

    %disp('-------------------------------------------');
    %disp(['Plotting for DUALCOLORINDEX=' num2str(DUALCOLORINDEX) ', TOPLOTFIELDNAMES=' TOPLOTFIELDNAMES{DUALCOLORINDEX}{plotFieldIndex}]);
    
    % Set plotfield
    TOPLOTFIELDNAME=TOPLOTFIELDNAMES{DUALCOLORINDEX}{plotFieldIndex};
    % Set dataset
    DATASETSTOPLOT = UBERLISTDATASETSTOPLOT{UberListIndex}; 
    
    % Run plotting script
    plottingGeneralDynamicData
    
end
end
end

%% ========================================================================
% Part II: What follows are plots for dual color experiments 
% =========================================================================

%% Plots per dataset, per parameter of interest ===========================
% Setting up

SELECTIONFIELD = 'groupID';
DATASETSTOPLOT = {'chromoCRP_cAMP800'};
TOPLOTFIELDNAMES = {{'concentrationCorrData', 'rateCorrData','concentrationDualCrossCorrData', 'rateDualCrossCorrData'},...
                    {'concentrationCorrData', 'rateCorrData'}}; 

ADDAVERAGETOFIRSTPLOT=1;      
CLEARFIGURES = 0;

if exist('myfig1','var'), if ~ishandle(myfig1), myfig1=figure(); else figure(myfig1); end, end
if exist('myfig2','var'), if ~ishandle(myfig2), myfig2=figure(); else figure(myfig2); end, end
if exist('myfig3','var'), if ~ishandle(myfig3), myfig3=figure(); else figure(myfig3); end, end
if exist('myfig4','var'), if ~ishandle(myfig4), myfig4=figure(); else figure(myfig4); end, end
                
%% Plotting of cross-correlation and related functions

for DUALCOLORINDEX=1:2
for plotFieldIndex=1:numel(TOPLOTFIELDNAMES{DUALCOLORINDEX})

    disp('-------------------------------------------');
    disp(['Plotting for DUALCOLORINDEX=' num2str(DUALCOLORINDEX) ', TOPLOTFIELDNAMES=' TOPLOTFIELDNAMES{DUALCOLORINDEX}{plotFieldIndex}]);        
    
    % Set plotfield
    TOPLOTFIELDNAME=TOPLOTFIELDNAMES{DUALCOLORINDEX}{plotFieldIndex};
        
    % Run plotting script
    plottingGeneralDynamicData
    
end
end

%% ========================================================================
% PART III: Kiviet style plots
% =========================================================================
%
% Plots per condition, with combined data for parameters of interest, that 
% is E-mu and P-mu plotted in the same graph, as in the Kiviet 2014 paper.
%
% Parameters are set separately per plot in each section below. There is
% some redundancy, which is also for clearity.
%
% Note that for the 2016-09-06 and the 2016-09-20 datasets, the 1st fluor
% is c, mCerulean (s70 reporter), and the 2nd fluor is y, mVenus (reporter
% for CRP activity). 
% This is conveyed to the script with the parameter DUALCOLORINDICES (e.g.
% setting DUALCOLORINDICES=[1,1] will result in two plots of Cerulean
% data).
%
%
%CONFIGFILE='U:\ZZ_EXPERIMENTAL_DATA\Data_per_project\CRPcAMP\config_projectCRPcAMP.m'

% Setting up
global crossCorrData;

SUBDIR = 'overview\';
SELECTIONFIELD = 'groupID';

if ~exist([PLOTOUTDIR SUBDIR],'dir')
    mkdir([PLOTOUTDIR SUBDIR]);
end

%% Loading
if ~exist('STOPRELOADING','var')
    plottingGeneral_sub_loaddataset
    STOPRELOADING = 1;
end

%% Set parameters

% Some need to be reset, because also used in previous sections
TOPLOTFIELDNAMES    = {'concentrationCorrData', 'rateCorrData'}; % this is always the case

if ~exist('GROUPSTOPLOT','var')
    GROUPSTOPLOT={'CRP_s70_chromosomal','chromoCRP_cAMP800','chromoCRP_cAMPLOW80','chromoCRP_cAMPHIGH5000'};
end
if ~exist('FILENAMESperGroupPartI','var')
    FILENAMESperGroupPartI = {'WT_','MED_','LOW_','HIGH_'};
end
if ~exist('FILENAMESperGroupPartII','var')
    FILENAMESperGroupPartII = {'CONST_','RCRP_'}; % Note that color1=const, color2 = rcrp
end

if ~exist('TITLESforgroupsPartI','var')
    TITLESforgroupsPartI = {'Wild type, ','\Delta cAMP, cAMP=800uM,  ','\Delta cAMP, cAMP=80uM,  ','\Delta cAMP, cAMP=5000uM, '};
end
if ~exist('TITLESforgroupsPartII','var')
    TITLESforgroupsPartII = {'const. reporter','CRP reporter'};
end

Line1Colors = linspecer(numel(GROUPSTOPLOT));
Line2Colors = ones(numel(GROUPSTOPLOT),3)*0; % also allows gray

%% Actually plotting the cross-corrs Daan style: do it using loop

optionsStruct=struct;
gatheredCCs = struct;
for groupIdx = 1:numel(GROUPSTOPLOT)
    for dualColorIdx = 1:2
        
        % parameters that are used by plotting script
        DATASETSTOPLOT = {GROUPSTOPLOT{groupIdx} GROUPSTOPLOT{groupIdx}};
        DUALCOLORINDICES =  [dualColorIdx dualColorIdx];
        LINECOLORS = {Line1Colors(groupIdx,:), Line2Colors(groupIdx,:)};

        % plotting script
        optionsStruct.STOPRELOADING=1;
        [hCC,output]=plottingGeneral_v2_CCs(DATASETSTOPLOT,DUALCOLORINDICES,LINECOLORS,SELECTIONFIELD,TOPLOTFIELDNAMES,optionsStruct); 

        % save average lines in struct also
        gatheredCCs.(GROUPSTOPLOT{groupIdx}).data = output;        
        
        % plot cosmetics
        figure(hCC.Number); 
        title([TITLESforgroupsPartI{groupIdx},TITLESforgroupsPartII{dualColorIdx}]);
        xlim([-10*output.hrsPerDoublingMean,10*output.hrsPerDoublingMean]);
        
        % Save file adjusted xlims
        for extIdx = 1:numel(EXTENSIONS)
            fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_xlimsadj_' FILENAMESperGroupPartI{groupIdx} FILENAMESperGroupPartII{dualColorIdx} '.' EXTENSIONS{extIdx}];
            if strcmp(EXTENSIONS{extIdx},'eps'), saveas(hCC,fileName,'epsc'); else saveas(hCC,fileName); end
        end
        
        xlim([-10,10]);
        
        % Save file with same xlims
        for extIdx = 1:numel(EXTENSIONS)
            fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_' FILENAMESperGroupPartI{groupIdx} FILENAMESperGroupPartII{dualColorIdx} '.' EXTENSIONS{extIdx}];
            if strcmp(EXTENSIONS{extIdx},'eps'), saveas(hCC,fileName,'epsc'); else saveas(hCC,fileName); end
        end
    end
end

disp('Done plotting Cross corrs');
winopen([PLOTOUTDIR SUBDIR]);


%% ========================================================================
% PART III-B. Now make make probability densities etc.
% =========================================================================

BINS=50;
GROUPSTOPLOT={'CRP_s70_chromosomal','chromoCRP_cAMP800','chromoCRP_cAMPLOW80','chromoCRP_cAMPHIGH5000'};
LEGENDNAMES ={'Feedback','Medium','Low','High'};
%mycolors=colorblind(2:end-1,:);
mycolors=linspecer(numel(GROUPSTOPLOT));

% Load the database w. datasets
if ~exist('DONTLOAD','var')
    %CONFIGFILE='U:\ZZ_EXPERIMENTAL_DATA\Data_per_project\CRPcAMP\config_projectCRPcAMP.m'
    plottingGeneral_sub_loaddataset
end

% But don't load it in the plottingGeneral_v2_PDF script. (Or next time
% called.)
DONTLOAD=1;

% What to plot (see plottingGeneral_v2_PDF for more info).
FIELDSOFINTEREST =...
        {{'muP9_fitNew_all','muP15_fitNew_all','muP23_fitNew_all','muP27_fitNew_all'},...
         {'C6_mean'}, ...
         {'Y6_mean'},...
         {'dC5'},...
         {'dY5'}};
TIMEFIELDS =...
        {'time',...
         'time_atC', ...
         'time_atY',...
         'dC5_time',...
         'dY5_time'};     
ALTERNATIVEFIELDNAMES = ...
        {'Growth rate',...
         '[Constitutive reporter]', ...
         '[CRP.cAMP reporter]',...
         'Production rate constitutive reporter',...
         'Production rate CRP.cAMP reporter'};         
SELECTIONFIELD='groupID';
SHOWRAWPLOTS=0;

% Output (sub)directory
SUBDIR = 'overview_PDF\';
if ~exist([PLOTOUTDIR SUBDIR],'dir'), mkdir([PLOTOUTDIR SUBDIR]), end % create if doesnt exist

disp('Section done');

%% Now gather info per dataset..
gatheredOutput = struct;
for idx=1:numel(GROUPSTOPLOT)
    % which one
    DATASETSTOPLOT = {GROUPSTOPLOT{idx}};
    % make pdf
    plottingGeneral_v2_PDF
    % store output
    gatheredOutput.(GROUPSTOPLOT{idx}) = output;
        % note this will have the structure gatheredOutput.(somefieldname).somedata{datasetIdx}{fieldIdx}
end
   

disp('Info gathering done');

%% Now plot the PDFs, and calculate the noises

% go over fields of interest
for paramOfInterestIdx = 1:numel(FIELDSOFINTEREST)

    figure; clf; hold on;

    legendLines=[];
    for groupIndex = 1:numel(GROUPSTOPLOT)

        currentGroup=GROUPSTOPLOT{groupIndex};
        
        % go over datasets
        for i = 1:numel(gatheredOutput.(currentGroup).allData)        

            %{
            % determine timewindow
            timeWindowSelectionIndices = ...
                [gatheredOutput.(currentGroup).allTimes{i}{paramOfInterestIdx} > gatheredOutput.(currentGroup).info{i}{paramOfInterestIdx}.fitTime(1)] & ...
                [gatheredOutput.(currentGroup).allTimes{i}{paramOfInterestIdx} < gatheredOutput.(currentGroup).info{i}{paramOfInterestIdx}.fitTime(2)];            
            %}
            
            % select current data
            theCurrentDataForTimeWindow = [gatheredOutput.(currentGroup).allDataForTimeWindow{i}{paramOfInterestIdx}];
            
            % get probability distribution function (PDF), ie hist
            [n,c]=hist(theCurrentDataForTimeWindow,BINS);
                        
            % normalize pdf
            dc=c(2)-c(1);
            area=sum(n)*dc;
            nNorm=n/area;            
            
            % plot
            l=plot(c,nNorm,'-','Color',mycolors(groupIndex,:),'LineWidth',2);
            
        end
       
        legendLines(end+1)=l;
        
    end
    
    % cosmetics & labels
    legend(legendLines,LEGENDNAMES);
    xlabel([info{1}{paramOfInterestIdx}.fieldType]);
    ylabel(['Normalized probability']);
    title([gatheredOutput.(currentGroup).info{1}{paramOfInterestIdx}.Fields{1} 10 ...
               ALTERNATIVEFIELDNAMES{paramOfInterestIdx}]...,
          ,'Interpreter','None');
    MW_makeplotlookbetter(15);

    % Save file
    for extIdx = 1:numel(EXTENSIONS)
        fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_pdf_' gatheredOutput.(currentGroup).info{1}{paramOfInterestIdx}.Fields{1} '.' EXTENSIONS{extIdx}];
        if strcmp(EXTENSIONS{extIdx},'eps'), saveas(gcf,fileName,'epsc'); else saveas(gcf,fileName); end
    end
    
end

%% Plot noises

% go over fields of interest
barvalues={};
for paramOfInterestIdx = 1:numel(FIELDSOFINTEREST)

    figure; clf; hold on;

    legendLines=[]; barcounter=0; barvalues{paramOfInterestIdx}=[];
    for groupIndex = 1:numel(GROUPSTOPLOT)

        currentGroup=GROUPSTOPLOT{groupIndex};
        
        % go over datasets
        for i = 1:numel(gatheredOutput.(currentGroup).allData)        
            
            barcounter=barcounter+1;
            
            % plot
            theNoiseValue=gatheredOutput.(currentGroup).noises{i}{paramOfInterestIdx};
            l=bar(barcounter,theNoiseValue,'FaceColor',mycolors(groupIndex,:));%,'LineWidth',2);
            barvalues{paramOfInterestIdx}(barcounter)=theNoiseValue;
            
        end
       
        legendLines(end+1)=l;
        
    end
    
    % cosmetics & labels
    ylim([0, max(barvalues{paramOfInterestIdx})*1.2]);
    legend(legendLines,LEGENDNAMES);
    xlabel('Replicates');
    ylabel(['Cv (std. dev./mean)']);
    title([gatheredOutput.(currentGroup).info{1}{paramOfInterestIdx}.Fields{1} 10 ...
               ALTERNATIVEFIELDNAMES{paramOfInterestIdx}]...,
          ,'Interpreter','None');
    MW_makeplotlookbetter(15);
    
    % Save file
    for extIdx = 1:numel(EXTENSIONS)
        fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_noise_' gatheredOutput.(currentGroup).info{1}{paramOfInterestIdx}.Fields{1} '.' EXTENSIONS{extIdx}];
        if strcmp(EXTENSIONS{extIdx},'eps'), saveas(gcf,fileName,'epsc'); else saveas(gcf,fileName); end
    end
    
end

%% ========================================================================
% Now make make branch plot overview etc.
% =========================================================================

% (Still using params from part III-B.)
% > GROUPSTOPLOT, mycolors, LEGENDNAMES

SUBDIR = 'overview_branches\';
if ~exist([PLOTOUTDIR SUBDIR],'dir')
    mkdir([PLOTOUTDIR SUBDIR]);
end

whichParamsToPlot = ...
            {{'X_time',    'muP9_fitNew_cycCor'},...
             {'dX5_time',  'dX5_cycCor'},...
             {'X_time',    'X6_mean_cycCor'}};

theYLabels = {'Growth [dbls/hr]',...
         'Production fluor X [a.u.]'...
         'X Fluor concentration [a.u.]'};
     
whichParamsToPlot2  = {{'dX5_time',  'dX5_cycCor'},{'X_time',    'X6_mean_cycCor'}};
theYLabels2         = {'Production fluor X [a.u.]','X Fluor concentration [a.u.]'};
        
winopen([PLOTOUTDIR SUBDIR]);

%% Start calculating and plotting

for groupIdx=1:numel(GROUPSTOPLOT)

    %% find datasets w. this identifier
    dataIdxs = find(strcmp({crossCorrData(:).groupID},GROUPSTOPLOT{groupIdx}));

    for dataIdx=dataIdxs

        %% Load schnitz file
        disp('Loading schnitz file');
        load(crossCorrData(dataIdx).schnitzfile);

        %%
        p.fitTime = crossCorrData(dataIdx).params.fitTime;

        % fluor1
        h1=MW_plotting_branches_mystyle(p,schnitzcells,mycolors(groupIdx,:),whichParamsToPlot,'Y',theYLabels);
        % fluor2
        h2=MW_plotting_branches_mystyle(p,schnitzcells,mycolors(groupIdx,:),whichParamsToPlot2,'C',theYLabels2);

        figure(h1); 
        supertitle([LEGENDNAMES{groupIdx} ' (id=' num2str(dataIdx) ')'],'FontSize',15,'Color',mycolors(groupIdx,:));
        figure(h2); 
        supertitle([LEGENDNAMES{groupIdx} ' (id=' num2str(dataIdx) ')'],'FontSize',15,'Color',mycolors(groupIdx,:));
        
        %crossCorrData(dataIdx).params
    
        % Save files
        for extIdx = 1:numel(EXTENSIONS)
            % fluor 1
            fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_branch_' LEGENDNAMES{groupIdx} '_Y_' num2str(dataIdx) '_' GROUPSTOPLOT{groupIdx} '.' EXTENSIONS{extIdx}];
            if strcmp(EXTENSIONS{extIdx},'eps'), saveas(h1,fileName,'epsc'); else saveas(h1,fileName); end
            % fluor 2
            fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_branch_' LEGENDNAMES{groupIdx} '_C_'  num2str(dataIdx) '_' GROUPSTOPLOT{groupIdx} '.' EXTENSIONS{extIdx}];
            if strcmp(EXTENSIONS{extIdx},'eps'), saveas(h2,fileName,'epsc'); else saveas(h2,fileName); end
        end
        
    end        
    
end

%% Now do the same thing again but make scatter plots..
% -------------------------------------------------------------------------

SUBDIR = 'overview_scatters\';
if ~exist([PLOTOUTDIR SUBDIR],'dir')
    mkdir([PLOTOUTDIR SUBDIR]);
end

NRCONTOURLINES=2;
FONTSIZE=15;

setsToScatterPlot = ...
    {   {'Y6_mean_cycCor','muP9_fitNew_cycCor'}, ...
        {'C6_mean_cycCor','muP9_fitNew_cycCor'}, ...
        {'dY5_cycCor','muP9_fitNew_atdY5_cycCor'}, ...
        {'dC5_cycCor','muP9_fitNew_atdC5_cycCor'}};

   
setLabels = ...    
    {'CRP reporter concentration (a.u.)',...
     'Const. reporter concentration (a.u.)',...
     'CRP reporter production rate (a.u.)',...
     'Const. reporter production rate (a.u.)'}
    
%%    

%% Initialize figures
myFigs = struct;    
for setIdx = 1:numel(setsToScatterPlot)
    myFigs(setIdx).hScatter=figure(); clf; hold on;    
    myFigs(setIdx).contourLineHandles=[]; myFigs(setIdx).scatterGroupsFirstHandles=[];
    myFigs(setIdx).xmaxes=[];
end

% Load schnitzcells data and plot scatter data for each
count=0; totalToLoad=sum(ismember({crossCorrData(:).groupID},GROUPSTOPLOT));
theScatterMeans={};
for groupIdx=1:numel(GROUPSTOPLOT)

    %% find datasets w. this identifier
    dataIdxs = find(strcmp({crossCorrData(:).groupID},GROUPSTOPLOT{groupIdx}));

    for dataIdx=dataIdxs

        count=count+1;
        
        %% Load schnitz file
        disp(['Loading schnitz file (' num2str(count) '/' num2str(totalToLoad) ')']);
        load(crossCorrData(dataIdx).schnitzfile);       
        
        for setIdx = 1:numel(setsToScatterPlot)
    
            figure(myFigs(setIdx).hScatter.Number);
            
            param1=setsToScatterPlot{setIdx}{1};
            param2=setsToScatterPlot{setIdx}{2};

            %%
            datax=[schnitzcells.(param1)];
            datay=[schnitzcells.(param2)];

            myFigs(setIdx).xmaxes(end+1)=max(datax);
            
            % contour (from kde)
            [bandwidth,density,X,Y] = kde2d([datax;datay]');      
            [C, lC] = contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);
            % organize means
            % Some means are already stored in crossCorrData, but this is
            % more straightforward
            
            theScatterMeans{groupIdx}{dataIdx}{setIdx} = [mean(datax(~isnan(datax))), mean(datay(~isnan(datay)))];
            
            myFigs(setIdx).contourLineHandles(end+1)=lC;

            % For later storage
            %{
            bandwidths{delayIdx}    = bandwidth;
            densities{delayIdx}     = density;
            Xs{delayIdx}            = X;
            Ys{delayIdx}            = Y;
            %}        
        
            %%
            myFigs(setIdx).lS = plot(datax,datay,'.','Color',mycolors(groupIdx,:));
            
        end                    
        
    end        
    
    for setIdx = 1:numel(setsToScatterPlot)
        myFigs(setIdx).scatterGroupsFirstHandles(end+1)=myFigs(setIdx).lS;
    end
    
end

disp('section done');

%% Now go over the figures once again and adjust cosmetics, then save
for setIdx = 1:numel(setsToScatterPlot)
    
    % Go to figure 
    currentH=figure(myFigs(setIdx).hScatter.Number);
    
    % Set contourlines on top
    uistack(myFigs(setIdx).contourLineHandles, 'top')
    
    % Plot mean values on top
    for groupIdx=1:numel(GROUPSTOPLOT)
        dataIdxs = find(strcmp({crossCorrData(:).groupID},GROUPSTOPLOT{groupIdx}));
        for dataIdx=dataIdxs
            
            l=plot(theScatterMeans{groupIdx}{dataIdx}{setIdx}(1),theScatterMeans{groupIdx}{dataIdx}{setIdx}(2),'o');
            set(l,'LineWidth',3,'MarkerFaceColor',mycolors(groupIdx,:),'MarkerEdgeColor','k','MarkerSize',15);
            
        end
    end
    
    % plot legend
    [h,icons,plots,legend_text] = legend(myFigs(setIdx).scatterGroupsFirstHandles,LEGENDNAMES,'FontSize',FONTSIZE)
    
    % adjust legend icon size
    for k = 1:numel(icons)%6:2:12%numel(icons)
        try
            icons(k).MarkerSize = 30;
        end   
        try
            icons(k).FontSize = FONTSIZE; % is again necessary for some reason
        end 
    end
    
    % More cosmetics
    MW_makeplotlookbetter(FONTSIZE);
    xlabel(setLabels(setIdx));
    ylabel('Growth rate (dbl/hr)');

    xlim([0,max(myFigs(setIdx).xmaxes)]);
    ylim([0,2]);
    
    % and save
    % Save files
    for extIdx = 1:numel(EXTENSIONS)
        % fluor 1
        fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_scatter_' setsToScatterPlot{setIdx}{1} '.' EXTENSIONS{extIdx}];
        if strcmp(EXTENSIONS{extIdx},'eps'), saveas(currentH,fileName,'epsc'); else saveas(currentH,fileName); end        
    end
end

winopen([PLOTOUTDIR SUBDIR]);

%%

winopen(PLOTOUTDIR);
disp('One script to rule them all done.');      




