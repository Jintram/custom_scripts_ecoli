% Calculates to which part mu-related fluctuations are responsible for
% extrinsic and intrinsic noise
% Total extrinsic noise resp. intrinsic noise is split up according to law
% of total covariance resp. variance
% extr noise: Cov(C,Y|mu)= Cov(<C|mu>,<Y|mu>) + <Cov(C,Y|mu>) = explained + unexplained
% intr noise: Var(C-Y)   = Var(<C-Y|mu>)      + <Var(C-Y|mu)> = explained + unexplained

% ************************************************************
% This function is very similar to "ExtrNoise_ExplainedByMu.m" but there
% are two differences:
% 1) The input is not a schnitzcells file but the already extracted Vectors
%    (i.e. also already taken care of fitTime and rm schnitzes)
%    for growth & production rate (in [...]/Extrinsic_Noise/Data_Collection/Data/)
% 2) Bins can be chosen equidistant or with equal amounts of datapoints in
%    each bin
% ************************************************************
 
% ***************
% Main Ouput: 'ordereredData' (contains all information)
% ***************

% ******************
% NO HISTORY DEPENDENCE IMPLEMENTED. FOR ONE TIME STEP IN HISTORY, 
% use ExtrIntrNoise_ExplainedByMu_1TimeptHist.m
% ******************



% ************************************************************************
%% 0) FACULTATIVE: LOAD A .MAT FILE CONTAINING SCHNITZCELLS AND EXTRACTED VECTORS
% But schnitzname has to be set!
% ************************************************************************
% 
clear myschnitzname
 %myschnitzname='Maltose full 2012-05-08 pos2';  % Adjust
 myschnitzname='Maltose full 2012-04-24 pos5';  % Adjust
% load(['\\biofysicasrv\Users2\Walker\ExtrinsicNoise\Data_Collection\Data\' myschnitzname '.mat']);
%


% ************************************************************************
%% 1) SET PLOT/BINNING OPTIONS AND GET DATA TOGETHER
% ************************************************************************

% ----------------------------------------------------------
% ADJUST
% ----------------------------------------------------------
% Plotting options
PLOTMUHIST=0;   % mu histogram
PLOTCONDVALUESEXTR=0;   % extr noise: <C|mu> vs mu plots (same for Y and Cov)
PLOTCONDVALUESINTR=0;   % intr noise: <Y-C|mu> vs mu plots (same for Var)
PLOTBINDEPENDENCE=0;    % explained fraction vs # bins used
PLOTCONDVARCOV=0;       % BLUBB TO IMPROVE, see at end of file
%close all % (de)activate for new/overwriting figures

% Data used [careful: in old version class was 'char' now it's 'double']
%field1=dC5_cycCor;  % production rate 1. will be normalized to mean=1
%field1str='dC5_cycCor'; % as string % BLUBB Automate
%field2=dY5_cycCor;  % production rate 2. will be normalized to mean=1
%field2str='dY5_cycCor'; % as string
%field3=mu_Third_cycCor; % growth rate. all vectors need to have same length!
%field3str='mu_Third_cycCor'; % as string


% *********BLUBB ************
field1=CFPrnd';
field2=YFPrnd';
field3=muRnd';
field1str='CFPrnd';
field2str='YFPrnd';
field3str='Murnd';
myschnitzname='0.2 corr. cfp and yfp ca 0.035 corr via mu';
% *********BLUBB ************

% Normalization
NormalizeStdDev=0; % ** standard=0. if=1: field1 and field2 are normalized to mean=1,
                   %        stddev is not normalized
                   % ** if =1, normalizes the standard deviations of field1 and field2 to =1. 
                   %        sets mean to =0 (maybe way of normalization still needs to be improved!!)
                   %        Then: ALSO field3 (MU) is standardized!! (change?)
                   %        All fields get mean=0!!! (careful when binning)
                   % Currently only implemented for MU_BINS_EQUALDIST=1!!!

% Binning
MUBINS_EQUALDIST=0; % =1: MuBins are of equal size (equally spaced) (but different # datapoints)
                    %     Bin_center is the MEAN of mu_min & mu_max of a
                    %     bin (-> equally spaced)
                    % =0: MuBins have each the same #datapoints
                    %     Bin_center is the MEDIAN of all datapoints within
                    %     that bin
                    % [if any other value than =0: behaves as =1]
                    
% EQ: only relevant for MUBINS_EQUALDIST=1 (otherwise these variables are cleared)
% BLUBB IS THE MEASURE x*MEAN GOOD? it's more robust for outliers...
SpacingMuBins=[0.08 0.1 0.15];%[0.07 0.085 0.11 0.14 0.22 0.45];%1/20];  % size of one bin, in fraction of <mu>. can be array.
SpecialMuBinCondValues=[0.22]; %extract all conditional value arrays (e.g. <yfp|mu>) for specific mu-bin
                % into extra variables for better handling. This bin size has to be included in 'SpacingMuBins'
% Note: Correctness of binning has been tested with: only 2 bins resp. very
% wide spacing: 0% explained. One bin for every datapoint:  100%
% explained

% NOTE: only relevant for MUBINS_EQUALDIST=0 (otherwise these variables are cleared)
%               NumberMuBins ExplFrac: Is only correct if none of the bins
%               is empty (usually the case)
%NumberMuBins=[2:10 12:2:20 25:5:40 50:10:100 ceil(length(dY5_cycCor)/6) ceil(length(dY5_cycCor)/5) ...
%    ceil(length(dY5_cycCor)/4) ceil(length(dY5_cycCor)/3) ceil(length(dY5_cycCor)/2) ceil(length(dY5_cycCor)/1)];%[5 10 15 20 25 30];    % # of bins for mu. centered around median of the bin.
                    % can be vector
NumberMuBins=ceil(length(field1)./[1:200 210:10:2000]);    % # of bins for mu. centered around median of the bin.
SpecialNumberBinsCondValues=1600; %extract all conditional value arrays (e.g. <yfp|mu>) for specific # mu-bins
                % into extra variables for better handling. This # bins has to be included in 'NumberMuBins'
% *************************************************************************

% ----------------------------------------------------------
% PREPARE DATA
% ----------------------------------------------------------

% % create strings from input variables (to label figures)
%varname=@(x) inputname(1); %BLUBB
%field1str=varname(field1);
%field2str=varname(field2);
%field3str=varname(field3);

% determine which type of mu-binning is performed
if MUBINS_EQUALDIST==0
    clear SpacingMuBins SpecialMuBinCondValues % for security
    disp('........')
    disp('BINNING MODE: EQUAL NUMBER OF DATAPOINTS IN EACH BIN')
    disp('USED # BINS:')
    disp(NumberMuBins)
    disp('........')
else
    MUBINS_EQUALDIST=1; % set to =1 (if not happened yet)
    clear NumberMuBins SpecialNumberBinsCondValues
    disp('........')
    disp('BINNING MODE: EQUALLY SPACED BINS')
    disp('USED MU BIN SIZES (Fraction of <mu>):')
    disp(SpacingMuBins)
    disp('........')
end

% combine data (a relict from the old version based on schnitzcells)
mymatrix=[field1',field2',field3']; %i.e.: [prodrate1, prodrate2, growthrate]

%normalize fluorescence fields (field1 & field2)
mymatrix(:,1)=mymatrix(:,1)./mean(mymatrix(:,1));
mymatrix(:,2)=mymatrix(:,2)./mean(mymatrix(:,2));
% normalize stddev if option is set (-> perform zscore)
if NormalizeStdDev
    if MUBINS_EQUALDIST==1
        mymatrix(:,1)=zscore(mymatrix(:,1));
        mymatrix(:,2)=zscore(mymatrix(:,2));
        % get mu bin sizes before subtracting mean. divide by std deviation
        % (->zscore rescaling)
        MuBinsAbsStd=(SpacingMuBins.*mean(mymatrix(:,3)))./std(mymatrix(:,3));
        mymatrix(:,3)=zscore(mymatrix(:,3));
    else
        error('Error: zscore normalization not implemented for non-equally spaced mu-bins')
    end
    
end

% % if some growth rates are <0 set artificially to=0 (should hardly ever happen)
% % necessary for consistent mu binning. only for non-extra normalization necessary.
% % alternative solution: delete this data.
%if ~NormalizeStdDev
%    idx=find(mymatrix(:,3)<0);
%    mymatrix(mymatrix(:,3)<0,3)=0;
%    disp(['# datapoints with mu<0 that were artificially set to=0: ' num2str(length(idx))]);
%end

%
% % ********** possibility to restrict size of data points (mymatrix) ***
%cut1=00001;
%cut2=60000;
%cut11=max(cut1,1); cut22=min(cut2,size(mymatrix,1));
%if cut1>1 | cut22<size(mymatrix,1)
%    disp('  ')
%    disp('------------------------------------------------------------------')
%    disp('Careful!! You artificially restricted to a subset of datapoints! (Set ''cut1'',''cut2'')')
%end
%mymatrix=mymatrix(cut11:cut22,:);
%clear cut1 cut2 cut11 cut22
% % ***************************

%get average mu
mumean=mean(mymatrix(:,3));

clear idx varname



% ************************************************************************
%% 2) CREATE STRUCT THAT CONTAINS DATA BINNED ACCORDING TO MU (orderedData)
% ************************************************************************

% create structure array 'orderedData' with length as number of elements in 
%    'SpacingMuBins' resp. 'NumberMuBins' (equally spaced vs equally populated).
% orderedData(i) contains all analysis for bin spacing SpacingMuBins(i)
%   (e.g. 0.2) resp. NumberBins(i) (e.g. 10) [both now called "binspacing (i)"]
orderedData=struct([]);
if MUBINS_EQUALDIST==1
    for i=1:length(SpacingMuBins)
        orderedData(i).SpacingMuBins_Current=SpacingMuBins(i);
        orderedData(i).Data=struct([]); % Data (specified further below)
        orderedData(i).NumberMuBins_Current=[];
        orderedData(i).type='equally spaced';
    end
else
    for i=1:length(NumberMuBins)
        orderedData(i).NumberMuBins_Current=NumberMuBins(i);
        orderedData(i).Data=struct([]); % Data (specified further below)
        orderedData(i).SpacingMuBins_Current=NaN;
        orderedData(i).type='equally populated';
    end
end

% -----------------------------------------------------------
% loop over all binning sizes (i) and fill the orderedData(i).Data structure
% -----------------------------------------------------------
% tempData:=orderedData(i).Data     -> all results for binspacing {i}
% tempData(k)  -> all results for k'th bin for a specific binspacing(i))
%                that is:
%            - muMin: lower border for mu for which data points land in this
%                   container. 
%                   for equally spaced: = (k-1)*(i)*<mu>    (i*<mu> is the size of one bin)
%            - muMax: upper border of this bins. 
%                   for equally spaced: = (k)*(i)*<mu>
%            - muCenter: centered mu of this bin. 
%                   for equally spaced: In the middle of muMin and muMax (mean)
%                   for equally populated: Median of all mu datapoints
%                          within this bin
%            - subMatrixMuBin: submatrix of mymatrix for which mu is within the borders
%            - probMuBin: fraction of all data points that lie within this muBin

% -----------------------------------------------------
% PERFORM THE MU BINNING
% if not standardized variables (especially mu!)
% -----------------------------------------------------
if ~NormalizeStdDev
    % ---------------------------
    % Loop over all binspacings (i)
    % ---------------------------
    for i=1:size(orderedData,2) % same size as SpacingMuBins or NumberMuBins
        tempData=struct([]);
        % ---------------------------
        % BINS EQUALLY SPACED
        % ---------------------------
        if MUBINS_EQUALDIST==1
            % determine absolute size of one Bin (in dbl/h)
            BinSizeAbs=SpacingMuBins(i)*mumean;
            % determine necessary number of bins for the given spacing
            % (bins used to start at mu=0, is changed now 2014-12)
             NumBins=ceil((max(mymatrix(:,3))-min(mymatrix(:,3)))/BinSizeAbs);
             % maybe Round is Better than ceil? but should only matter in
             % "artefact" cases with very many bins (what are disadvantages
             % of round?)
            orderedData(i).NumberMuBins_Current=NumBins;
            % global minimum of mu (beginning of first bin)
            muMinGlobal=min(mymatrix(:,3));
            % loop over each mu bin and extract data that lies within this bin
            for n=1:NumBins
                tempData(n).muMin=(n-1)*BinSizeAbs+muMinGlobal; 
                tempData(n).muMax=n*BinSizeAbs+muMinGlobal;
                tempData(n).muCenter=(n-0.5)*BinSizeAbs+muMinGlobal;
                idx=find(mymatrix(:,3)>=tempData(n).muMin & mymatrix(:,3)<tempData(n).muMax);
                tempData(n).subMatrixMuBin=mymatrix(idx,:);
                tempData(n).probMuBin=length(idx)./size(mymatrix,1);
            end
        else
        % ---------------------------
        % BINS EQUALLY POPULATED
        % ---------------------------
            % determine # datapoints per Bin
            NumDatapoints=size(mymatrix,1);
            NumDatapointsBin=ceil(NumDatapoints/NumberMuBins(i));
            % sort mymatrix according to growth rate 
            [mymatrixsorted,sortindex]=sortrows(mymatrix,3);
            % loop over each mu bin and extract data that lies within this bin
            for n=1:NumberMuBins(i)
                % idx of first and last datapoint
                idxmin=min(1+(n-1)*NumDatapointsBin,NumDatapoints);
                idxmax=min(n*NumDatapointsBin,NumDatapoints);
                % growth rates (mumin&mumax: if not border bins: take
                %    average of last point in this bin and first point in
                %    previous/next bin
                switch n
                    case 1
                        tempData(n).muMin=mymatrixsorted(1,3);
                        idxnext=idxmax+1;
                        if idxnext>size(mymatrixsorted,1), idxnext=idxmax; end
                        tempData(n).muMax=0.5*(mymatrixsorted(idxmax,3)+mymatrixsorted(idxnext,3));
                        tempData(n+1).muMin=tempData(n).muMax;
                    case NumberMuBins(i) | NumDatapoints % not sure if | works
                        tempData(n).muMax=mymatrixsorted(end,3);
                    otherwise
                        idxnext=idxmax+1;
                        if idxnext>size(mymatrixsorted,1), idxnext=idxmax; end
                        tempData(n).muMax=0.5*(mymatrixsorted(idxmax,3)+mymatrixsorted(idxnext,3));
                        if n<NumberMuBins(i) % this is a random statement
                            tempData(n+1).muMin=tempData(n).muMax;
                        end
                end
                % sorted or not
                idx=find(mymatrix(:,3)>=tempData(n).muMin & mymatrix(:,3)<tempData(n).muMax+0.00000001); % the 0....01 makes sure the highest datapoint is not ignored
                tempData(n).subMatrixMuBin=mymatrix(idx,:);
                tempData(n).muCenter=median(mymatrix(idx,3));                
                tempData(n).probMuBin=length(idx)./size(mymatrix,1);
            end
        end
        orderedData(i).Data=tempData;
        
        % --------------------
        % Plot histogram of growth rates
        % --------------------
        if PLOTMUHIST
            % get data together
            xdata=zeros(length(tempData),1);
            ydata=zeros(length(tempData),1);
            for run=1:length(xdata)
                xdata(run)=tempData(run).muCenter;
                ydata(run)=tempData(run).probMuBin;
            end
            figure
            set(gcf,'WindowStyle','docked')
            clf
            bar(xdata,ydata)
            %hold on
            title('growth rate distribution');
            xlabel('mu')
            ylabel('fraction')
            if MUBINS_EQUALDIST==1
                    myannotation=['Bins eq spaced. ' num2str(SpacingMuBins(i)) '*<mu>'];
            else
                    myannotation=['Bins eq populated. ' num2str(NumberMuBins(i)) ' bins.'];
            end
            annotation(gcf,'textbox',...
             [0.661714285714286 0.680952380952381 0.193642857142857 0.145238095238098],...
             'String',{['<mu>=' num2str(mumean)],myannotation}, 'FitBoxToText','on');
        end

    end % loop over all binspacings. no extra normalization
    
% -----------------------------------------------------
% PERFORM THE MU BINNING
% if extra standardized distributions
% -----------------------------------------------------
else
    if MUBINS_EQUALDIST==1
        
        for i=1:length(SpacingMuBins)
            tempData=struct([]);
            % absolute size of one Bin (in dbl/h)
            BinSizeAbs=MuBinsAbsStd(i);
            % get necessary number of bins to capture also highest mu. (bins start
            % at mu<0
            NumBins=ceil(max(mymatrix(:,3))/BinSizeAbs-min(mymatrix(:,3))/BinSizeAbs);
            tempData(n).NumberMuBins_Current=NumBins;
            MuStart=min(mymatrix(:,3));
            % loop over each mu bin and extract data that lies within this bin
            for n=1:NumBins
                tempData(n).muMin=MuStart+(n-1)*BinSizeAbs;
                tempData(n).muMax=MuStart+n*BinSizeAbs;
                tempData(n).muCenter=MuStart+(n-0.5)*BinSizeAbs;
                idx=find(mymatrix(:,3)>=tempData(n).muMin & mymatrix(:,3)<tempData(n).muMax);
                tempData(n).subMatrixMuBin=mymatrix(idx,:);
                tempData(n).probMuBin=length(idx)./size(mymatrix,1);
            end
            orderedData(i).Data=tempData;

            % --------------------
            % Plot histogram of growth rates
            % --------------------
            if PLOTMUHIST
                % get data together
                xdata=zeros(length(tempData),1);
                ydata=zeros(length(tempData),1);
                for run=1:length(xdata)
                    xdata(run)=tempData(run).muCenter;
                    ydata(run)=tempData(run).probMuBin;
                end
                figure
                set(gcf,'WindowStyle','docked')
                clf
                bar(xdata,ydata)
                %hold on
                title(['growth rate distribution']);
                xlabel('mu')
                ylabel('fraction')
                if MUBINS_EQUALDIST==1
                    myannotation=['Bins eq spaced. ' num2str(SpacingMuBins(i)) '*<mu>'];
                else
                    myannotation=['Bins eq populated. ' num2str(NumberMuBins(i)) ' bins.'];
                end
                annotation(gcf,'textbox',...
                 [0.661714285714286 0.680952380952381 0.193642857142857 0.145238095238098],...
                 'String',{['<mu>=' num2str(mumean)],myannotation}, 'FitBoxToText','on');
            end

        end % loop over all binspacings. with extra normalization   
    else
        error('Error: zscore normalization not implemented for non-equally spaced mu-bins')
    end
    
end % if ~NormalizeStdDev

clear BinSizeAbs NumBins i idx n run tempData xdata ydata myannotation muMinGlobal
        

% ************************************************************************
%% 3) CALCULATE CONDITIONAL (CO)VARIANCES AND EXPECTATION VALUES IN EACH BIN
% ************************************************************************
% i.e.: calculate them separately for each mu-Bin (and of course separately
% for different binning sizes)

%loop over all binspacings (i)
for i=1:size(orderedData,2)  
    % for readability:
    tempData=orderedData(i).Data;
    % loop over each mu bin
    for n=1:length(tempData)
        %for readability:
        f1=tempData(n).subMatrixMuBin(:,1); %1st field (e.g. cfp production rate)
        f2=tempData(n).subMatrixMuBin(:,2);
        % all rows of subMatrixMuBin are equally likely
        
        % --------------------------
        % Extrinsic: conditional covariance Cov(C,Y|mu)
        %            conditional expectations <C|mu>, <Y,mu>
        % --------------------------
        tempData(n).CondCov_field1field2=mean(f1.*f2)-mean(f1)*mean(f2); %normalized by N (cf help cov)
        %if n==10  %debug
        %    mean(f1.*f2)-mean(f1)*mean(f2)
        %    cov(f1,f2,1)
        %end
        tempData(n).CondMean_field1=mean(f1);
        tempData(n).CondMean_field2=mean(f2);
        
        % --------------------------
        % Intrinsic: conditional variance Var(C-Y|mu)
        %            conditional expectation <C-Y|mu>
        % --------------------------
        Diff_f1f2=f1-f2;  % more accurate name: Diff_field1field2, but would be very long
        tempData(n).Diff_f1f2=Diff_f1f2;
        tempData(n).CondMean_Diff_f1f2=mean(Diff_f1f2);
        tempData(n).CondVar_Diff_f1f2=mean(Diff_f1f2.*Diff_f1f2)-mean(Diff_f1f2)^2; %normalized by N (cf help cov)
        
    end
    orderedData(i).Data=tempData;
    
    % ------------------------------------
    % Plot conditional expectations/covariance for extr noise
    % ------------------------------------
    if PLOTCONDVALUESEXTR
        % get data together 
        xdata=zeros(length(tempData),1);
        covdata=zeros(size(xdata));
        meanf1data=zeros(size(xdata));
        meanf2data=zeros(size(xdata));
        
        for run=1:length(xdata)
            xdata(run)=tempData(run).muCenter;
            covdata(run)=tempData(run).CondCov_field1field2;
            meanf1data(run)=tempData(run).CondMean_field1;
            meanf2data(run)=tempData(run).CondMean_field2;
        end
        figure
        set(gcf,'WindowStyle','docked')
        clf
        plot(xdata,covdata,'.-r')
        hold on
        plot(xdata,meanf1data,'.-b')
        plot(xdata,meanf2data,'.-g')
        % plot the edges of the bins
        muBinEdges=unique([[orderedData(i).Data.muMin],[orderedData(i).Data.muMax]]);
        plot(muBinEdges,zeros(size(muBinEdges)),'.k','MarkerSize',10)
        
        if MUBINS_EQUALDIST==1
           mytitle=['Extrinsic: Conditional Means & Cov. Bins eq. spaced ' num2str(SpacingMuBins(i)) '*<mu>.'];
        else
            mytitle=['Extrinsic: Conditional Means & Cov. Bins eq. populated ' num2str(NumberMuBins(i)) ' bins.'];
        end
        title(mytitle);
        xlabel('mu (muCenter)')
        ylabel('<..|mu> resp. cov(..|mu)')
        legend('Cov(field1,field2|mu)', '<field1|mu>', '<field2|mu>','Location','NW')
        grid on
    end
    
    % ------------------------------------
    % Plot conditional expectations/variance for intr noise
    % ------------------------------------
    if PLOTCONDVALUESINTR
        % get data together 
        xdata=zeros(length(tempData),1);
        vardata=zeros(size(xdata));
        meanDiff_f1f2data=zeros(size(xdata));
        
        for run=1:length(xdata)
            xdata(run)=tempData(run).muCenter;
            vardata(run)=tempData(run).CondVar_Diff_f1f2;
            meanDiff_f1f2data(run)=tempData(run).CondMean_Diff_f1f2;
        end
        figure
        set(gcf,'WindowStyle','docked')
        clf
        plot(xdata,vardata,'.-r')
        hold on
        plot(xdata,meanDiff_f1f2data,'.-b')
        % plot the edges of the bins
        muBinEdges=unique([[orderedData(i).Data.muMin],[orderedData(i).Data.muMax]]);
        plot(muBinEdges,zeros(size(muBinEdges)),'.k','MarkerSize',10)
        
        if MUBINS_EQUALDIST==1
           mytitle=['Intrinsic: Conditional Means & Var. Bins eq. spaced ' num2str(SpacingMuBins(i)) '*<mu>.'];
        else
            mytitle=['Intrinsic: Conditional Means & Var. Bins eq. populated ' num2str(NumberMuBins(i)) ' bins.'];
        end
        title(mytitle)
        xlabel('mu (muCenter)')
        ylabel('<...|mu> resp. var(..|mu)')
        legend('Var((field1-field2)|mu)', '<(field1-field2)|mu>','Location','NW')
        grid on
        
    end
    
end % loop over all binspacings

clear Diff_f1f2    f1 f2 i meanDiff_f1f2data 
%clear meanf1data meanf2data n run tempData vardata xdata



% ************************************************************************
%% 4) Calculate explained & unexplained parts in extr & intr noise
% ************************************************************************
% law of total (co)variance

% loop over binspacings
for i=1:size(orderedData,2)  
    % --------------------
    % create vectors with conditional means,cov,var of all bins
    % --------------------
    
    % for readability:
    tempData=orderedData(i).Data;
    % get data together
    % -------------------
    % extr noise
    % -------------------
    muCenter_vec=[];%zeros(length(tempData),1);  % in plots named: xdata
    CondCov_field1field2_vec=[];%zeros(size(muCenter));
    CondMean_field1_vec=[];%zeros(size(muCenter));
    CondMean_field2_vec=[];%zeros(size(muCenter));
    probMuBin_vec=[];%zeros(size(muCenter));
    % -------------------
    % intr noise
    % -------------------
    CondVar_Diff_f1f2_vec=[];
    CondMean_Diff_f1f2_vec=[];
    
    % loop over each bin and extract data. ONLY USE NOT-NAN data (if data
    % contains NaN, then this bin is empty and has probability =0 -> will
    % not contribute to averages in an analytical solution)
    for n=1:length(tempData);
        if ~isnan(tempData(n).CondCov_field1field2) % bin is not empty
            %muCenter_vec necessary?
            CondCov_field1field2_vec(end+1)=tempData(n).CondCov_field1field2;
            CondMean_field1_vec(end+1)=tempData(n).CondMean_field1;
            CondMean_field2_vec(end+1)=tempData(n).CondMean_field2;
            probMuBin_vec(end+1)=tempData(n).probMuBin;
            CondVar_Diff_f1f2_vec(end+1)=tempData(n).CondVar_Diff_f1f2;
            CondMean_Diff_f1f2_vec(end+1)=tempData(n).CondMean_Diff_f1f2;
        end
    end
    
    % --------------------
    % Extrinsic Noise: calculate law of total covariance
    % --------------------
    Mean_of_Cov_Unexplained_Extr=sum(probMuBin_vec.*CondCov_field1field2_vec);
    % better readability
    x=CondMean_field1_vec; y=CondMean_field2_vec;
    Cov_of_Mean_Explained_Extr=sum((x.*y).*probMuBin_vec)-sum(x.*probMuBin_vec)*sum(y.*probMuBin_vec);
    %  Cov(<x|mu>,<y|mu>) = Av_over_mu(<x|mu>*<y|mu>) -  Av_over_mu(<x|mu>) * Av_over_mu(<y|mu>)
    
    orderedData(i).Mean_of_Cov_Unexplained_Extr=Mean_of_Cov_Unexplained_Extr;
    orderedData(i).Cov_of_Mean_Explained_Extr=Cov_of_Mean_Explained_Extr;
    
    % add total covariance - independent of binning, because calculated
    % directly:
    orderedData(i).totalCov_Extr=mean(mymatrix(:,1).*mymatrix(:,2))-mean(mymatrix(:,1))*mean(mymatrix(:,2));
    
    % ******
   % % add alternative way of calculation BLUBB -> gave absolutely same
   % result
   % orderedData(i).Mean_of_Cov_Unexpl_Extr_ALT=mean(mymatrix(:,1).*mymatrix(:,2))-sum((x.*y).*probMuBin_vec);
   % orderedData(i).Cov_of_Mean_Expl_Extr_ALT=sum((x.*y).*probMuBin_vec)-mean(mymatrix(:,1))*mean(mymatrix(:,2));
    % ******
    
    % --------------------
    % Intrinsic Noise: calculate law of total variance
    % --------------------
    % better readability
    v=CondVar_Diff_f1f2_vec; m=CondMean_Diff_f1f2_vec;
    Mean_of_Var_Unexplained_Intr=sum(probMuBin_vec.*v);
    Var_of_Mean_Explained_Intr=sum((m.*m).*probMuBin_vec)-sum(m.*probMuBin_vec)^2;
    
    orderedData(i).Mean_of_Var_Unexplained_Intr=Mean_of_Var_Unexplained_Intr;
    orderedData(i).Var_of_Mean_Explained_Intr=Var_of_Mean_Explained_Intr;
    
    % add total variance of (f1-f2) - independent of binning, because calculated
    % directly:
    diff=mymatrix(:,1)-mymatrix(:,2);
    orderedData(i).totalVar_Intr=mean(diff.*diff)-mean(diff)^2;
       
    % ******
    % add alternative way of calculation BLUBB  -> gave absolutely same
    % result
    %orderedData(i).Mean_of_Var_Unexpl_Intr_ALT=orderedData(i).totalVar_Intr-sum(((x-y).^2).*probMuBin_vec);
    %orderedData(i).Var_of_Mean_Expl_Intr_ALT=sum(((x-y).^2).*probMuBin_vec)-(mean(mymatrix(:,1)-mymatrix(:,2)))^2;
    % ******
    
end % loop over binspacings

% ------------------------------------------
% Report the results. The actual calculation is finished here
% ------------------------------------------

disp(' ')
disp('_______________________________________________________________')
disp(['Analyzed ' myschnitzname '.'])
disp('_______________________________________________________________')
disp(['Inputs were extra normalized (zscore); 1=true,0=false : ' num2str(NormalizeStdDev)]);
disp(['Used fields ' field1str ' & ' field2str ])
disp(['and conditioned on ' field3str ' without history.'])
disp(['Spacing of Bins: 1=equally spaced, 0=equally populated:   ' num2str(MUBINS_EQUALDIST)])
disp('  ')
disp('Extrinsic')
if MUBINS_EQUALDIST==1
    disp(['BinFrac(of <mu>)    #bins (calc)     Explained_Extr  Unexplained_Extr   Total        Frac_Expl     Frac_Unexpl'])
    for i=1:length(orderedData)
        disp([num2str(SpacingMuBins(i)) '                      ' num2str(orderedData(i).NumberMuBins_Current) '             ' ...
            num2str(orderedData(i).Cov_of_Mean_Explained_Extr) '         ' ... 
            num2str(orderedData(i).Mean_of_Cov_Unexplained_Extr) ...
            '         '   num2str(orderedData(i).totalCov_Extr)  '        '  ...
            num2str(orderedData(i).Cov_of_Mean_Explained_Extr/orderedData(i).totalCov_Extr) ...
            '       ' num2str(orderedData(i).Mean_of_Cov_Unexplained_Extr/orderedData(i).totalCov_Extr)]);
    end
else
    disp(['BinFrac(of <mu>)    #bins            Explained_Extr  Unexplained_Extr   Total        Frac_Expl     Frac_Unexpl'])
    for i=1:length(orderedData)
        disp(['variable             ' num2str(orderedData(i).NumberMuBins_Current) '                ' ...
            num2str(orderedData(i).Cov_of_Mean_Explained_Extr) '        ' ... 
            num2str(orderedData(i).Mean_of_Cov_Unexplained_Extr) ...
            '            '   num2str(orderedData(i).totalCov_Extr)  '      '  ...
            num2str(orderedData(i).Cov_of_Mean_Explained_Extr/orderedData(i).totalCov_Extr) ...
            '       ' num2str(orderedData(i).Mean_of_Cov_Unexplained_Extr/orderedData(i).totalCov_Extr)]);
    end
end
% *******Alternative*********** BLUBB
%disp(['BinFrac     Explained_Extr  Unexplained_Extr   Total        Frac_Expl     Frac_Unexpl ALTERNATIVE'])
%for i=1:length(orderedData)
%    disp([num2str(SpacingMuBins(i)) '         ' num2str(orderedData(i).Cov_of_Mean_Expl_Extr_ALT) '        ' ... 
%        num2str(orderedData(i).Mean_of_Cov_Unexpl_Extr_ALT) ...
%        '            '   num2str(orderedData(i).totalCov_Extr)  '      '  ...
%        num2str(orderedData(i).Cov_of_Mean_Expl_Extr_ALT/orderedData(i).totalCov_Extr) ...
%        '       ' num2str(orderedData(i).Mean_of_Cov_Unexpl_Extr_ALT/orderedData(i).totalCov_Extr)]);
%end
% ******************

disp(' ')
disp('Intrinsic')
if MUBINS_EQUALDIST
    disp(['BinFrac(of <mu>)    #bins (calc)     Explained_Intr  Unexplained_Intr   Total        Frac_Expl     Frac_Unexpl'])
    for i=1:length(orderedData)
        disp([num2str(SpacingMuBins(i)) '                    ' num2str(orderedData(i).NumberMuBins_Current) '                ' ...
            num2str(orderedData(i).Var_of_Mean_Explained_Intr) '      ' ... 
            num2str(orderedData(i).Mean_of_Var_Unexplained_Intr) ...
            '         '   num2str(orderedData(i).totalVar_Intr)  '     '  ...
            num2str(orderedData(i).Var_of_Mean_Explained_Intr/orderedData(i).totalVar_Intr) ...
            '       ' num2str(orderedData(i).Mean_of_Var_Unexplained_Intr/orderedData(i).totalVar_Intr)]);
    end
else
    disp(['BinFrac(of <mu>)    #bins            Explained_Intr  Unexplained_Intr   Total        Frac_Expl     Frac_Unexpl'])
    for i=1:length(orderedData)
        disp(['variable             ' num2str(orderedData(i).NumberMuBins_Current) '                 ' ...
            num2str(orderedData(i).Var_of_Mean_Explained_Intr) '      ' ... 
            num2str(orderedData(i).Mean_of_Var_Unexplained_Intr) ...
            '           '   num2str(orderedData(i).totalVar_Intr)  '     '  ...
            num2str(orderedData(i).Var_of_Mean_Explained_Intr/orderedData(i).totalVar_Intr) ...
            '       ' num2str(orderedData(i).Mean_of_Var_Unexplained_Intr/orderedData(i).totalVar_Intr)]);
    end
end

% *******Alternative*********** BLUBB -> lead to exactly same result
% (implementation thus correct...)
%disp(['BinFrac     Explained_Intr  Unexplained_Intr   Total        Frac_Expl     Frac_Unexpl ATLTERNATIVE'])
%for i=1:length(orderedData)
%    disp([num2str(SpacingMuBins(i)) '         ' num2str(orderedData(i).Var_of_Mean_Expl_Intr_ALT) '        ' ... 
%        num2str(orderedData(i).Mean_of_Var_Unexpl_Intr_ALT) ...
%        '           '   num2str(orderedData(i).totalVar_Intr)  '     '  ...
%        num2str(orderedData(i).Var_of_Mean_Expl_Intr_ALT/orderedData(i).totalVar_Intr) ...
%        '       ' num2str(orderedData(i).Mean_of_Var_Unexpl_Intr_ALT/orderedData(i).totalVar_Intr)]);
%end
% ******************

disp('_______________________________________________________________')

clear CondCov_field1field2_vec CondMean_Diff_f1f2_vec CondMean_field1_vec
clear CondMean_field2_vec CondVar_Diff_f1f2_vec diff i n muCenter_vec m probMuBin_vec
clear tempData v x y



% *******************************************************************
%% 4) INVESTIGATE DEPENDENCE OF EXPLAINED FRACTION ON #BINS
% *******************************************************************

% create matrix with: fracbinSize --- #bins --- fracexplExtr --- fracexplIntr
matBinDependence=zeros(size(orderedData,2),4);
if MUBINS_EQUALDIST==1
    matBinDependence(:,1)=SpacingMuBins;
else
    matBinDependence(:,1)=NaN;
end
matBinDependence(:,2)=[orderedData.NumberMuBins_Current];
matBinDependence(:,3)=[orderedData.Cov_of_Mean_Explained_Extr]./[orderedData.totalCov_Extr];
matBinDependence(:,4)=[(orderedData.Var_of_Mean_Explained_Intr)]./[(orderedData.totalVar_Intr)];

if MUBINS_EQUALDIST==0
    disp(' ')
    disp('Some plots may stay empty because size of bins is not fix for MUBIN_EQUALDIST=0.')
    disp(' ')
end
if PLOTBINDEPENDENCE
    figure('Position',[100 100 800 800])
    clf
    hold on
    subplot(2,2,1)
    plot(matBinDependence(:,1),matBinDependence(:,3),'.-' , 'LineWidth',2,'MarkerSize',15)
    xlabel('fracBinSize')
    ylabel('frac_explained','Interpreter','None')
    title('Extrinsic Noise','FontSize',12)
    yy=get(gca,'ylim');
    if yy(2)>0  % for 'normal' data always the case
        set(gca,'ylim',[0 yy(2)]);
    end
    grid on
    if ~NormalizeStdDev
        subplot(2,2,3)
        plot(matBinDependence(:,2),matBinDependence(:,3),'.-' , 'LineWidth',2,'MarkerSize',15)
        %xlim([0 40])
        xlabel('#bins')
        ylabel('frac_explained','Interpreter','None');
        set(gca, 'XDir','reverse')
        if yy(2)>0  % for 'normal' data always the case
          set(gca,'ylim',[0 yy(2)]);
       end
    grid on
    end
    subplot(2,2,2)
    plot(matBinDependence(:,1),matBinDependence(:,4),'.-' , 'LineWidth',2,'MarkerSize',15)
    xlabel('fracBinSize')
    ylabel('frac_explained','Interpreter','None')
    title('Intrinsic Noise','FontSize',12)
    yy=get(gca,'ylim'); 
    if yy(2)>0  % for 'normal' data always the case
        set(gca,'ylim',[0 yy(2)]);
    end
    grid on
    if ~NormalizeStdDev
        subplot(2,2,4)
        plot(matBinDependence(:,2),matBinDependence(:,4),'.-' , 'LineWidth',2,'MarkerSize',15)
        %xlim([0 40])
        xlabel('#bins')
        ylabel('frac_explained','Interpreter','None');
        set(gca, 'XDir','reverse')
        if yy(2)>0  % for 'normal' data always the case
           set(gca,'ylim',[0 yy(2)]);
        end
        grid on
    end
    annotation(gcf,'textbox',[0.191 0.485 0.82525 0.0412500000000001],...
        'String',{['Plotted ' num2str(size(mymatrix,1)) ' Schnitzes of ' myschnitzname '.'], ...
        [field1str ', ' field2str ' conditioned on ' field3str '.']},...
        'FitBoxToText','on','Interpreter','none', 'FontSize',12);
end

clear  meanf1data meanf2data run xdata vardata covdata

%********** BLUBB
figure
clf; hold on;
set(gcf,'WindowStyle','docked')
plot(matBinDependence(:,2),matBinDependence(:,3),'.-' , 'LineWidth',2,'MarkerSize',15)
%xlim([0 40])
xlabel('#bins')
ylabel('frac_explained','Interpreter','None');
%set(gca, 'XDir','reverse')
yy=get(gca,'ylim');
if yy(2)>0  % for 'normal' data always the case
     set(gca,'ylim',[0 yy(2)]);
end
grid on
title([myschnitzname '  ' field1str  ' ' field2str  '  ' field3str '  . extr noise'])




% **********************************************************************
%% 5) Extract conditional averages and variances for specific mu-bin into new
% variables   
% **********************************************************************

CondMean_f1=[]; % conditional mean of field 1 (e.g. yfprate)
CondMean_f2=[];
CondCov_f1f2=[];
CondMean_f1minf2=[];   % conditional mean of field1 - field2   :   <yfp-cfp|mu>
CondVar_f1minf2=[];    % conditional Variance of field1 - field2   :   Var(yfp-cfp|mu)
MuVec=[];    % vector which contains the mean of each mu-bin
ProbMuBin=[];

% find the binspacing of interest
if MUBINS_EQUALDIST==1
    idxbinspacing=find([orderedData.SpacingMuBins_Current]==SpecialMuBinCondValues);
    if isempty(idxbinspacing)
        disp(' ')
        disp('Special Spacing is not in Calculated Dataset (evaluating subcell (5)). Will return.')
        return
    end
else
    idxbinspacing=find([orderedData.NumberMuBins_Current]==SpecialNumberBinsCondValues);
    if isempty(idxbinspacing)
        disp(' ')
        disp('Special Spacing is not in Calculated Dataset (evaluating subcell (5)). Will return.')
        return
        return
    end
end
i=idxbinspacing;
subData=orderedData(1,i).Data;
for runmu=1:length(subData) % loop over each mu-bin
            CondMean_f1=[CondMean_f1, subData(runmu).CondMean_field1]; % conditional mean of field 1 (e.g. yfprate)
            CondMean_f2=[CondMean_f2, subData(runmu).CondMean_field2];
            CondCov_f1f2=[CondCov_f1f2, subData(runmu).CondCov_field1field2];
            CondMean_f1minf2=[CondMean_f1minf2, subData(runmu).CondMean_Diff_f1f2];   % conditional mean of field1 - field2   :   <yfp-cfp|mu>
            CondVar_f1minf2=[CondVar_f1minf2, subData(runmu).CondVar_Diff_f1f2];
            MuVec=[MuVec, subData(runmu).muCenter];
            ProbMuBin=[ProbMuBin, subData(runmu).probMuBin];
end
% check for empty bins (-> NaNvalues)
idx=find(ProbMuBin>0);
CondMean_f1=CondMean_f1(idx);
CondMean_f2=CondMean_f2(idx);
CondCov_f1f2=CondCov_f1f2(idx);
CondMean_f1minf2=CondMean_f1minf2(idx);
CondVar_f1minf2=CondVar_f1minf2(idx);
MuVec=MuVec(idx);

%BLUBB
if PLOTCONDVARCOV
    % DO YOU WANT TO NORM??
    figure
    set(gcf,'WindowStyle','docked')
    clf
    % ConVar_Norm:  normalized by <yfp|mu>*<cfp|mu>
    CondVar_Norm=CondVar_f1minf2./(CondMean_f1.*CondMean_f2);
    hold on
    %CondVar_Norm=CondVar_Norm/mean(CondVar_Norm); 
    
    %plot(MuVec,CondVar_f1minf2/mean(CondVar_f1minf2),'.-','Color',[0 0.8 1])
    plot(MuVec,CondVar_f1minf2,'.-','Color',[0 0.8 1])
    plot(MuVec,CondVar_Norm,'.-m')
    legend('CondVar (same as red line in prev. plot)','CondVar normed w average of bin: <f1|mu><f2|mu>')
    %title('Conditional Variance (Var(f1-f2|mu). Intr noise. Normed to mean=1')
    title('TOIMPROVE. Conditional Variance (Var(f1-f2|mu). Intr noise. no norm')
    xlabel(['mu (binfrac:' num2str(SpecialMuBinCondValues)])
    ylabel('var(...|mu)')
    grid on
end

clear i runmu idx PLOTCONDVARCOV