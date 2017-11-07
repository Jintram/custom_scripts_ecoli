
%% Parameters for all code ================================================

% CONFIGFILE='U:\ZZ_EXPERIMENTAL_DATA\Data_per_project\CRPcAMP\config_projectCRPcAMP.m';

warning('TODO: move the files in the output folder (also .docx files) and change path here..');
PLOTOUTDIR = '\\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\A_Step5_Plots\CRPcAMP\';
EXTENSIONS = {'png','svg','fig'};

some_colors;

% =========================================================================

%%
% Measure timing
tic;

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
TOPLOTFIELDNAMES = {'concentrationCorrData', 'rateCorrData'}; % this is always the case

if ~exist('GROUPSTOPLOT','var')
    GROUPSTOPLOT={'CRP_s70_chromosomal','chromoCRP_cAMP800','chromoCRP_cAMPLOW80','chromoCRP_cAMPHIGH5000'};
end

if ~exist('FLUORCOLORS','var')
    FLUORCOLORS={'C','Y'};
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
SUBDIR = 'overview_ccs\';
if ~exist([PLOTOUTDIR SUBDIR],'dir')
    mkdir([PLOTOUTDIR SUBDIR]);
end

optionsStruct=struct;
gatheredCCs = struct;
for groupIdx = 1:numel(GROUPSTOPLOT)
    for dualColorIdx = 1:2
        
        % The following code plots multiple cross-correlations into one
        % plot.
        
        % parameters that are used by plotting script
        % ===
        % the datasets the CCs come from; here usually the same
        currentDataSetsToPlot = {GROUPSTOPLOT{groupIdx} GROUPSTOPLOT{groupIdx}};
        % # of the two parameters that the cross-corr is calculated for
        DUALCOLORINDICES =  [dualColorIdx dualColorIdx];
        % colors of the two lines
        LINECOLORS = {Line1Colors(groupIdx,:), Line2Colors(groupIdx,:)};

        % plotting script
        optionsStruct.STOPRELOADING=1;
        [hCC,output]=plottingGeneral_v2_CCs(currentDataSetsToPlot,DUALCOLORINDICES,LINECOLORS,SELECTIONFIELD,TOPLOTFIELDNAMES,optionsStruct); 

        % save average lines in struct also
        gatheredCCs.(GROUPSTOPLOT{groupIdx}).data.(FLUORCOLORS{dualColorIdx}) = output;        
        
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

%% Repeat for autocorrelations
TOPLOTFIELDNAMESAUTOCORR={'concentrationautoCorrData', 'rateautoCorrData','growthautoData'};
TITLESAUTOCORR={'Concentration','Rate','Growth rate'};
LEGENDNAMES={'E-E','p-p','\mu-\mu'};
% note that handling of autocorrelations for growth rates is a bit awkward
% as it is done twice for fluor1 and fluor2

optionsStructAC=struct;
gatheredACs = struct;
for groupIdx = 1:numel(GROUPSTOPLOT)
    for dualColorIdx = 1:2
        
        % parameters that are used by plotting script
        currentDataSetsToPlot = {GROUPSTOPLOT{groupIdx} GROUPSTOPLOT{groupIdx} GROUPSTOPLOT{groupIdx}};
        DUALCOLORINDICES =  [dualColorIdx dualColorIdx dualColorIdx];
        LINECOLORS = {Line1Colors(groupIdx,:), Line2Colors(groupIdx,:) [.5 .5 .5]};

        % plotting script
        optionsStructAC.STOPRELOADING=1;
        optionsStructAC.LEGENDNAMES=LEGENDNAMES;
        [hAC,output]=plottingGeneral_v2_CCs(currentDataSetsToPlot,DUALCOLORINDICES,LINECOLORS,SELECTIONFIELD,TOPLOTFIELDNAMESAUTOCORR,optionsStructAC); 

        % save average lines in struct also
        gatheredACs.(GROUPSTOPLOT{groupIdx}).data.(FLUORCOLORS{dualColorIdx}) = output;        
        
        % plot cosmetics
        figure(hAC.Number); 
        title([TITLESforgroupsPartI{groupIdx},TITLESforgroupsPartII{dualColorIdx}]);
        xlim([0,10*output.hrsPerDoublingMean]);              
        
    end
end

%% Repeat for CCs between two colors 
TOPLOTFIELDNAMESCOLORCORR={'concentrationDualCrossCorrData', 'rateDualCrossCorrData'};
TITLESAUTOCORR={'Concentration','Rate'};
LEGENDNAMES={'E_{consti.}-E_{CRP}','p_{consti.}-p_{CRP}'};
% note that handling of autocorrelations for growth rates is a bit awkward
% as it is done twice for fluor1 and fluor2

fieldNamesCheckPerGroup={};
optionsStructCY=struct;
gatheredCYs = struct;
for groupIdx = 1:numel(GROUPSTOPLOT)
    for dualColorIdx = 1
        
        % parameters that are used by plotting script
        currentDataSetsToPlot = {GROUPSTOPLOT{groupIdx} GROUPSTOPLOT{groupIdx}};
        DUALCOLORINDICES =  [dualColorIdx dualColorIdx];
        LINECOLORS = {Line1Colors(groupIdx,:), Line2Colors(groupIdx,:)};

        % plotting script
        optionsStructCY.STOPRELOADING=1;
        optionsStructCY.LEGENDNAMES=LEGENDNAMES;
        [hCY,output]=plottingGeneral_v2_CCs(currentDataSetsToPlot,DUALCOLORINDICES,LINECOLORS,SELECTIONFIELD,TOPLOTFIELDNAMESCOLORCORR,optionsStructCY); 

        % save average lines in struct also
        gatheredCYs.(GROUPSTOPLOT{groupIdx}).data.CY = output;        
        
        % plot cosmetics
        figure(hCY.Number); 
        title([TITLESforgroupsPartI{groupIdx}]);
        xlim([-10*output.hrsPerDoublingMean,10*output.hrsPerDoublingMean]);                  
            
    end            
    
    % Tell user X-order of fields for R(X1,X2)
    fieldNamesCheckPerGroup{groupIdx} = crossCorrData(output.idxsCrossCorrData).rateDualCrossCorrFieldNames;
            
end

fieldNamesCheckPerGroup
disp('It might be a good idea to double-check the order of your fieldnames, since this depends on the order in the config. file.');
% so far it has been R(C,Y), so R(const., CRP)

%% Repeat for production-concentration correlation cross-correlations

TOPLOTFIELDNAMESCOLORCORR={'rateConcentrationCC'};
TITLESPECORR={'Wild type','Optimal cAMP (800uM)','Low cAMP (80uM)','High cAMP (5000uM)'};
TITLESPECORRp2 = {', Constitutive',', CRP cAMP'};
LEGENDNAMES={'Concentration-rate'};
% note that handling of autocorrelations for growth rates is a bit awkward
% as it is done twice for fluor1 and fluor2

fieldNamesCheckPerGroup={};
optionsStructPE=struct;
gatheredPEs = struct;
for groupIdx = 1:numel(GROUPSTOPLOT)
    for dualColorIdx = 1:2
        
        % parameters that are used by plotting script
        currentDataSetsToPlot = {GROUPSTOPLOT{groupIdx}};
        DUALCOLORINDICES =  [dualColorIdx];
        LINECOLORS = {Line1Colors(groupIdx,:)};

        % plotting script
        optionsStructPE.STOPRELOADING=1;
        optionsStructPE.LEGENDNAMES=LEGENDNAMES;
        [hCY,output]=plottingGeneral_v2_CCs(currentDataSetsToPlot,DUALCOLORINDICES,LINECOLORS,SELECTIONFIELD,TOPLOTFIELDNAMESCOLORCORR,optionsStructPE); 

        % save average lines in struct also
        gatheredPEs.(GROUPSTOPLOT{groupIdx}).data.(FLUORCOLORS{dualColorIdx}) = output;        
        
        % plot cosmetics
        figure(hCY.Number); 
        title([TITLESPECORR{groupIdx} TITLESPECORRp2{dualColorIdx}]);
        xlim([-10*output.hrsPerDoublingMean,10*output.hrsPerDoublingMean]);                  
            
    end            
    
    % Tell user X-order of fields for R(X1,X2)
    fieldNamesCheckPerGroup{groupIdx} = crossCorrData(output.idxsCrossCorrData).rateDualCrossCorrFieldNames;
            
end

fieldNamesCheckPerGroup
disp('It might be a good idea to double-check the order of your fieldnames, since this depends on the order in the config. file.');
% so far it has been R(C,Y), so R(const., CRP)

%% Plot autocorrelations in combined plot

SUBDIR = 'overview_autocorr\';
TYPENAMES={'Conc_','Rate_','Growth_'};

if ~exist([PLOTOUTDIR SUBDIR],'dir')
    mkdir([PLOTOUTDIR SUBDIR]);
end

for typeIndex=1:3
for dualColorIdx=1:2

    %%
    hCombinedAC = figure(); clf; hold on; 
    lines=[]; alltaus=[]; allRS =[];
    for groupIdx= 1:numel(GROUPSTOPLOT)

        % title
        title([TITLESAUTOCORR{typeIndex} ' ' TITLESforgroupsPartII{dualColorIdx}]);
        if typeIndex==3 % for growth rate TITLESforgroupsPartII is redundant
            title([TITLESAUTOCORR{typeIndex}]);
        end

        tau=gatheredACs.(GROUPSTOPLOT{groupIdx}).data.(FLUORCOLORS{dualColorIdx}).datatau{typeIndex};
        R=gatheredACs.(GROUPSTOPLOT{groupIdx}).data.(FLUORCOLORS{dualColorIdx}).datacorrelation{typeIndex};
        hrsPerDoublingMean=gatheredACs.(GROUPSTOPLOT{groupIdx}).data.(FLUORCOLORS{dualColorIdx}).hrsPerDoublingMean;
            % note lines 1-3 correspond to 
            % TOPLOTFIELDNAMESAUTOCORR={'concentrationautoCorrData', 'rateautoCorrData','growthautoData'}

        lineColor=Line1Colors(groupIdx,:);
        lines(end+1)=plot(tau,R,'-','LineWidth',2,'Color',lineColor);

        % plot the mean hrs per doubling as a bar
        plot([0,hrsPerDoublingMean],[-0.03*groupIdx,-.03*groupIdx],'-','LineWidth',4,'Color', lineColor);

        % gather data pile for x and ylims
        alltaus=[alltaus tau];
        allRs  =[allRS R -0.03*groupIdx-.02];
    end
        
    myXlim=[0,max(alltaus)];
    xlim(myXlim);

    myYlim=[min(allRs),1.01];
    ylim(myYlim);

    plot([myXlim],[0,0],'-k');

    legend(lines,LEGENDNAMES);

    xlabel('Time (hrs)');
    ylabel('Autocorrelation');
    MW_makeplotlookbetter(20);
    
    % Save plot
    for extIdx = 1:numel(EXTENSIONS)
        fileName = [PLOTOUTDIR SUBDIR EXTENSIONS{extIdx} '_autocorr_' TYPENAMES{typeIndex} FILENAMESperGroupPartII{dualColorIdx} '.' EXTENSIONS{extIdx}];
        if strcmp(EXTENSIONS{extIdx},'eps'), saveas(hCombinedAC,fileName,'epsc'); else saveas(hCombinedAC,fileName); end
    end
    
end
end

disp('Done plotting Auto corrs');
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

%% Get a better estimate of the variance, also important for modeling

% GROUPSTOPLOT was set earlier to determine which groups to plot
% Also, mycolors=linspecer(numel(GROUPSTOPLOT));

YFIELDSVARIANCE=    {'muP9_fitNew_all'   , 'dY5_sum'            ,'dC5_sum', 'C6_mean' ,'Y6_mean'};
NAMESVARIANCE  =    {'Variance in growth', 'Variance CRP prod.' ,'Variance consti. prod.','Var. in CRP label','Var. consti. label'};
INDICESSELECTION =  {0                   , 'indices_at_dC'       , 'indices_at_dC','indices_at_C','indices_at_C'};
    % 0=no selection

gatheredVarianceOutput = struct; countVar=0;
for groupIdx=1:numel(GROUPSTOPLOT)

    %% find datasets w. this identifier
    currentGroup = GROUPSTOPLOT{groupIdx};
    dataIdxs = find(strcmp({crossCorrData(:).groupID},currentGroup));

    disp(['Looking at : ' currentGroup]);
    
    % Loading applicable schnitzes and continue:
    for dataSetsInGroupIdx = 1:numel(dataIdxs)

        currentDataIdx = dataIdxs(dataSetsInGroupIdx);
        
        % Load the schnitz file
        clear schnitzcells
        load(crossCorrData(currentDataIdx).schnitzfile,'schnitzcells');

        % Fix it if necessary with an additional field
        if ~isfield(schnitzcells, 'indices_at_C')
            disp('Editing saved schnitzcells struct.');
            schnitzcells = MW_addToSchnitzes_indices_atXdX(crossCorrData(currentDataIdx).schnitzfile, 'C');
        end

        % administration
        countVar=countVar+1;
        
        % Find the variance for different fields in the schnitzcells struct
        for fieldIndex = 1:numel(YFIELDSVARIANCE)    

            %%
            fieldOfInterest = YFIELDSVARIANCE{fieldIndex};

            allFrames    = [schnitzcells.frame_nrs];
            allTimes     = [schnitzcells.time];

            % This is necessary for fields that are not present for all
            % timepoints
            if ischar(INDICESSELECTION{fieldIndex})
                indicesToSelect= find([schnitzcells.(INDICESSELECTION{fieldIndex})]);
            else
                indicesToSelect= 1:numel(allFrames);
            end

            % get applicable frames    
            theFrames    = allFrames(indicesToSelect);
            uniqueFrames = sort(unique(theFrames));

            % get applicable times    
            theTimes     = allTimes(indicesToSelect);
            uniqueTimes  = sort(unique(theTimes));

            % Load the field of interest
            theY      = [schnitzcells.(fieldOfInterest)];

            % calculate variances per frame and store them
            count=0;
            variances = NaN(1,numel(uniqueFrames)); variancesNormalized = NaN(1,numel(uniqueFrames));
            for frameNr = uniqueFrames

                count=count+1;
                selectionIdxs = theFrames==frameNr;    
                growthThisFrame = theY(selectionIdxs);
                growthThisFrame = growthThisFrame(~isnan(growthThisFrame));
                growthThisFrameNormalized = growthThisFrame./mean(growthThisFrame);
                variances(count) = var(growthThisFrame);
                variancesNormalized(count) = var(growthThisFrameNormalized);
            end

            % calculate an estimated variance based on last 10 frames
            meanVariance  = mean(variancesNormalized(end-9:end));
            medianVariance= median(variancesNormalized(end-9:end));            
            stdVariance   = std(variancesNormalized(end-9:end));                                    

            h=figure(); clf; hold on;
            %bar(uniqueTimes/60,variancesNormalized);%,'-o','LineWidth',2)
            plot(uniqueTimes/60,variancesNormalized,'.','LineWidth',2);          
            title(['Population based variance' 10 'Median last 10 frms = ' num2str(medianVariance) ' +/- ' num2str(stdVariance) '.']);
            xlabel('Time (hrs)'); 
            ylabel([NAMESVARIANCE{fieldIndex}]);
            MW_makeplotlookbetter(20);

            
            gatheredVarianceOutput.(fieldOfInterest).meanVariances(countVar)       = ...
                meanVariance;
            gatheredVarianceOutput.(fieldOfInterest).medianVariances(countVar)     = ...
                medianVariance;
            gatheredVarianceOutput.(fieldOfInterest).stdVariances(countVar)        = ...
                stdVariance;
            gatheredVarianceOutput.(fieldOfInterest).NormalizedVariancesTimeseries{countVar} = ...
                [uniqueTimes;variancesNormalized];
            gatheredVarianceOutput.(fieldOfInterest).groupNames{countVar}          = ...
                currentGroup;
            gatheredVarianceOutput.(fieldOfInterest).groupIdx(countVar)          = ...
                groupIdx;
            
                
            %{
            gatheredVarianceOutput(groupIdx).(fieldOfInterest).meanVariance{dataSetsInGroupIdx} = ...
                meanVariance;
            gatheredVarianceOutput(groupIdx).(fieldOfInterest).medianVariance{dataSetsInGroupIdx} = ...
                medianVariance;
            gatheredVarianceOutput(groupIdx).(fieldOfInterest).stdVariance{dataSetsInGroupIdx} = ...
                stdVariance;
            gatheredVarianceOutput(groupIdx).(fieldOfInterest).NormalizedVariancesTimeseries{dataSetsInGroupIdx} = ...
                [uniqueTimes;variancesNormalized];
            gatheredVarianceOutput(groupIdx).(fieldOfInterest).groupName= currentGroup;
            %}
            
        end        
        
    end
end

disp('Variance analysis done..');

%% Now make a plot to compare them
disp([gatheredVarianceOutput.(YFIELDSVARIANCE{1}).groupNames]);
for fieldIdx = 1:5

    figure; clf; hold on;   

    %{
    legendLines=[]; barcounter=0; barvalues{fieldIdx}=[];
    for groupIndex = 1:numel(GROUPSTOPLOT)

        currentGroup=GROUPSTOPLOT{groupIndex};

        % go over datasets
        for i = 1:numel(gatheredOutput.(currentGroup).allData)        

            barcounter=barcounter+1;

            % plot
            theNoiseValue=gatheredOutput.(currentGroup).noises{i}{fieldIdx};
            l=bar(barcounter,theNoiseValue,'FaceColor',mycolors(groupIndex,:));%,'LineWidth',2);
            barvalues{fieldIdx}(barcounter)=theNoiseValue;

        end

        legendLines(end+1)=l;

    end
    %}

    % plot
    legendLines=[]; lastGroupIdx=0;
    for i=1:numel(gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).medianVariances)

        X=i;
        Y=gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).medianVariances(i);
        %Y=gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).meanVariances(i);
        errY=gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).stdVariances(i);
        groupIdx=gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).groupIdx(i);
        theColor = mycolors(groupIdx,:);

        errorbar(X,Y,errY,'LineWidth',3,'Color','k');
        l=bar(X,Y,'FaceColor',theColor);    

        if (groupIdx~=lastGroupIdx)
            legendLines(end+1)=l;
        end
        lastGroupIdx=groupIdx;

    end

    % cosmetics
    barLocs    = 1:numel(gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).meanVariances);
    labelNames = gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).groupNames;
    set(gca,'TickLabelInterpreter','none'); 
    set(gca,'XTick',barLocs,'XTickLabel',labelNames);
    set(gca,'XTickLabelRotation',45)
    ylabel(NAMESVARIANCE{fieldIdx}); % Variance based on normalized growth
    MW_makeplotlookbetter(15);    

    legend(legendLines,LEGENDNAMES);
    xlabel('Replicates');

    disp([YFIELDSVARIANCE{fieldIdx} '='    num2str(gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).meanVariances)]);
    
end

disp('Section done');

%%
fieldIdx = 5;

legendLines=[]; lastGroupIdx=0;
figure; clf; hold on;
for i=1:numel(gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).medianVariances)
    groupIdx=gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).groupIdx(i);
    theColor = mycolors(groupIdx,:);
    dataLine = gatheredVarianceOutput.(YFIELDSVARIANCE{fieldIdx}).NormalizedVariancesTimeseries{i};
    l=plot(dataLine(1,:)./60,dataLine(2,:),'-','Color',theColor,'LineWidth',1)    
    
    if (groupIdx~=lastGroupIdx)
        legendLines(end+1)=l;
    end
    lastGroupIdx=groupIdx;
end

legend(legendLines,LEGENDNAMES);
xlabel('Time (hrs)');
ylabel(NAMESVARIANCE{fieldIdx});
MW_makeplotlookbetter(20);

ylim([0,0.7]);

%%

winopen(PLOTOUTDIR);
disp('One script to rule them all done.');      

t2=toc;
disp(['Done in ' num2str(t2/60) ' minutes']);




