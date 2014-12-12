% Calculates to which part mu-related fluctuations are responsible for
% extrinsic and intrinsic noise
% production rates (or concentrations) must be contained in field1 and
% field2 and will be normalized to average=1


% ***************
% Main Ouput: 'ordereredData' (contains all information)
% ***************

% ******************
% NO HISTORY DEPENDENCE IMPLEMENTED. FOR ONE TIME STEP IN HISTORY, 
% use ExtrIntrNoise_ExplainedByMu_1TimeptHist.m
% ******************

% **********************************************
% INITIALIZE AND EXTRACT DATA
% **********************************************

% *** FACULTATIVE ***
% run initiation and selection from schnitzcells excel file
%load '\\biofysicasrv\Users2\Walker\ExperimentsCollectedData\schnitzcellsCollectionForNoiseExplByMu20130510.mat';
% takes a lot of time!
%

%%

% **** ADJUST *****
PLOTMUHIST=1;
PLOTCONDVALUESEXTR=1;
PLOTCONDVALUESINTR=1;
PLOTBINDEPENDENCE=0;
PLOTCONDVARCOV=1;

% ** (de)activate **
close all
% ***********

schnitzUseName='schnitzcells_acetate20120620pos2_full';  %schnitzcells selection
%schnitzUseName='schnitzcells';
%schnitzUseName='schnitzcells2012_06_17_ASC631';  %schnitzcells selection
%field1='dY5_len_smooth3_dt_subtr';    % production rates, maybe concentrations
%field2='dC5_len_smooth3_dt_subtr';    % production rates, maybe concentrations
%field3='muP15_fitNew_subtr';         %growth rate!
%field1='dY5_len_dt_subtr';  % production rates, maybe concentrations
%field2='dC5_len_dt_subtr';  % production rates, maybe concentrations
%field3='muP15_fitNew_subtr'; %growth rate!
%field1='dY5_sum_dt_cycCor';
%field2='dC5_sum_dt_cycCor';
%field3='muP11_fitNew_cycCor';

field1='dC5_cycCor';  
field2='dY5_cycCor';
field3='muP15_fitNew_cycCor'; 

%field1='dR5_sum_dt_s_cycCor';  
%field2='dG5_sum_dt_s_cycCor';
%field3='muP15_fitNew_cycCor'; 
%field1='dR5_sum_dt_s_cycCor';  
%field2='dG5_sum_dt_s_cycCor';
%field3='muP15_fitNew_cycCor'; 


%field1='muP15_fitNew_cycCor';
%field2='muP15_fitNew_cycCor';
%field3='dY5_sum_dt_s_cycCor';
NormalizeStdDev=0; % ** if =1, normalizes the standard deviations of field1 and field2 to =1. 
                   %sets mean to =0 (maybe way of normalization still needs to be improved!!)
                   % Now: ALSO field3 (MU) is standardized!!
                   % All fields get mean=0!!! (careful when binning)
                   % ** if=0, field1 and field2 are normalized to mean=1,
                   % stddev is untouched.
%MuBins=[ 0.01 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6];%1/20];  % size of one bin, in fraction of <mu>. can be array.
MuBins=[0.2];%1/20];  % size of one bin, in fraction of <mu>. can be array.
SpecialMuBinCondValues=0.2; %extract all conditional value arrays (e.g. <yfp|mu>) for specific mu-bin
                % into extra variables for better handling. This bin size has to be included in 'MuBins'
% *****************

% get schitzcells structure (above: only name)
eval(['schnitzUse=' schnitzUseName ';']);
% extract data
mymatrix=zeros(0,3); % field1 - field2 - field3 (mu) (can include weights here)

% loop over all schnitzes
disp('Make sure you have the wanted useForPlot settings!')
for i=1:length(schnitzUse)
    s=schnitzUse(i);
    % use only schnitzes with useForPlot=1
    if ~existfield(s,'useForPlot') | s.useForPlot==1 %blubb
        
   % if s.useForPlot==1
        % check if fields actually contain data (length>0)
        % Note: e.g. prod.rates and mu do not always have same array length
        % (rates can be 1 shorter for late data)
        minLength=min([length(s.(field1)),length(s.(field2)),length(s.(field3))]);
        if minLength>0
            % check if NaN values exist
            Nan1=sum(isnan(s.(field1)(1:minLength)));
            Nan2=sum(isnan(s.(field2)(1:minLength)));
            Nan3=sum(isnan(s.(field3)(1:minLength)));
            if (Nan1+Nan2+Nan3)==0
                % add values to mymatrix
                addmatrix=[s.(field1)(1:minLength)',s.(field2)(1:minLength)',s.(field3)(1:minLength)'];
                mymatrix=[mymatrix;addmatrix];
            end
        end
    end
end
% ********** possibility to restrict size of data points (mymatrix) ***
cut1=00001;
cut2=60000;
cut11=max(cut1,1); cut22=min(cut2,size(mymatrix,1));
if cut1>1 | cut22<size(mymatrix,1)
    disp('  ')
    disp('------------------------------------------------------------------')
    disp('Careful!! You artificially restricted to a subset of datapoints! (Set ''cut1'',''cut2'')')
end
mymatrix=mymatrix(cut11:cut22,:);
clear cut1 cut2 cut11 cut22
% ***************************

%normalize
% fluocolors1 and 2
mymatrix(:,1)=mymatrix(:,1)./mean(mymatrix(:,1));
mymatrix(:,2)=mymatrix(:,2)./mean(mymatrix(:,2));
% zscore
if NormalizeStdDev
    mymatrix(:,1)=zscore(mymatrix(:,1));
    mymatrix(:,2)=zscore(mymatrix(:,2));
    % get mu bin sizes before subtracting mean. divide by std deviation
    % (->zscore rescaling)
    MuBinsAbsStd=(MuBins.*mean(mymatrix(:,3)))./std(mymatrix(:,3));
    mymatrix(:,3)=zscore(mymatrix(:,3));
    
end
% force negative growth rates to be =0 (blubb) if mu not extra normalized. Alternatively, this data
% could also be deleted. step necessary for consistent mu-binning (don't
% forget some data)
if ~NormalizeStdDev
    mymatrix(mymatrix(:,3)<0,3)=0;
end
 
clear Nan1 Nan2 Nan3 addmatrix i minLength s

%%
% **********************************************
% CREATE STRUCT THAT CONTAINS DATA BINNED ACCORDING TO MU (orderedData)
% **********************************************
% plotting options
%PLOTMUHIST=0;

%MuBins=[1/5,  1/10,    1/20,   1/100]; %redundant

%get average mu
mumean=mean(mymatrix(:,3));

%create structure array with length #MuBins that contains
% muBinSize (e.g. 0.2 or 0.1)) - in fraction of <mu>
% Data (specified further below)
orderedData=struct([]);
for i=1:length(MuBins)
    orderedData(i).muBinSize=MuBins(i);
    orderedData(i).Data=struct([]);
end

%loop over all binning sizes (i) and fill the orderedData(i).Data structure
% for readability:    tempData:=orderedData(i).Data
% tempData(k)  (k'th entry of tempData, and for a specific binning size
% (i)) contains all information about the k'th bin of mu if a binning size
% of (i) is used. That is
% - muMin: lower border for mu for which data points land in this
%          container. = (k-1)*(i)*<mu>    (i*<mu> is the size of one bin)
% - muMax: upper border of this bins. = (k)*(i)*<mu>
% - muAv: average (centered) mu of this bin. In the middle of muMin and
%         muMax
% - subMatrixMuBin: submatrix of mymatrix for which mu is within the borders
% - probMuBin: fraction of all data points that lie within this muBin

% -----------------------------------------------------
% if not standardized variables (especially mu!)
% -----------------------------------------------------
if ~NormalizeStdDev
    for i=1:length(MuBins)
        tempData=struct([]);
        % absolute size of one Bin (in dbl/h)
        BinSizeAbs=MuBins(i)*mumean;
        % get necessary number of bins to capture also highest mu. (bins start
        % at mu=0, all data with mu<0 is forced to mu=0)
        NumBins=ceil(max(mymatrix(:,3))/BinSizeAbs);
        % loop over each mu bin and extract data that lies within this bin
        for n=1:NumBins
            tempData(n).muMin=(n-1)*BinSizeAbs;
            tempData(n).muMax=n*BinSizeAbs;
            tempData(n).muAv=(n-0.5)*BinSizeAbs;
            idx=find(mymatrix(:,3)>=tempData(n).muMin & mymatrix(:,3)<tempData(n).muMax);
            tempData(n).subMatrixMuBin=mymatrix(idx,:);
            tempData(n).probMuBin=length(idx)./size(mymatrix,1);
        end
        orderedData(i).Data=tempData;

        if PLOTMUHIST
            % get data together
            xdata=zeros(length(tempData),1);
            ydata=zeros(length(tempData),1);
            for run=1:length(xdata)
                xdata(run)=tempData(run).muAv;
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
            annotation(gcf,'textbox',...
             [0.661714285714286 0.680952380952381 0.193642857142857 0.145238095238098],...
             'String',{['<mu>=' num2str(mumean)],['bin size frac=' num2str(MuBins(i))]}, 'FitBoxToText','on');
        end

    end % loop over all binSizes
    
% -----------------------------------------------------
% extra standardized distributions
% -----------------------------------------------------
else
    for i=1:length(MuBins)
        tempData=struct([]);
        % absolute size of one Bin (in dbl/h)
        BinSizeAbs=MuBinsAbsStd(i);
        % get necessary number of bins to capture also highest mu. (bins start
        % at mu<0
        NumBins=ceil(max(mymatrix(:,3))/BinSizeAbs-min(mymatrix(:,3))/BinSizeAbs);
        MuStart=min(mymatrix(:,3));
        % loop over each mu bin and extract data that lies within this bin
        for n=1:NumBins
            tempData(n).muMin=MuStart+(n-1)*BinSizeAbs;
            tempData(n).muMax=MuStart+n*BinSizeAbs;
            tempData(n).muAv=MuStart+(n-0.5)*BinSizeAbs;
            idx=find(mymatrix(:,3)>=tempData(n).muMin & mymatrix(:,3)<tempData(n).muMax);
            tempData(n).subMatrixMuBin=mymatrix(idx,:);
            tempData(n).probMuBin=length(idx)./size(mymatrix,1);
        end
        orderedData(i).Data=tempData;

        if PLOTMUHIST
            % get data together
            xdata=zeros(length(tempData),1);
            ydata=zeros(length(tempData),1);
            for run=1:length(xdata)
                xdata(run)=tempData(run).muAv;
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
            annotation(gcf,'textbox',...
             [0.661714285714286 0.680952380952381 0.193642857142857 0.145238095238098],...
             'String',{['<mu>=' num2str(mumean)],['bin size frac=' num2str(MuBins(i))]}, 'FitBoxToText','on');
        end

    end % loop over all binSizes
    
end     

clear BinSizeAbs PLOTMUHIST NumBins i idx n run tempData xdata ydata
        
        
 %%
% **********************************************
% CALCULATE CONDIITONAL (CO)VARIANCES AND EXPECTATION VALUES
% **********************************************   
% i.e.: calculate them seperately for each mu-Bin (and of course seperately
% for different binning sizes)
%PLOTCONDVALUESEXTR=0;
%PLOTCONDVALUESINTR=0;


% -----------
% conditional covariance
% ----------
%loop over all binning sizes (i)
for i=1:length(MuBins)  
    % for readability:
    tempData=orderedData(i).Data;
    % loop over each mu bin
    for n=1:length(tempData)
        %for readability:
        f1=tempData(n).subMatrixMuBin(:,1); %1st field (e.g. yfp production rate)
        f2=tempData(n).subMatrixMuBin(:,2);
        % all rows of subMatrixMuBin are equally likely
        
        % --- for extrinsic ---
        tempData(n).CondCov_field1field2=mean(f1.*f2)-mean(f1)*mean(f2); %normalized by N (cf help cov)
        %if n==10  %debug
        %    mean(f1.*f2)-mean(f1)*mean(f2)
        %    cov(f1,f2,1)
        %end
        tempData(n).CondMean_field1=mean(f1);
        tempData(n).CondMean_field2=mean(f2);
        
        % --- for intrinsic ---
        Diff_f1f2=f1-f2;  % more accurate name: Diff_field1field2, but would be very long
        tempData(n).Diff_f1f2=Diff_f1f2;
        tempData(n).CondMean_Diff_f1f2=mean(Diff_f1f2);
        tempData(n).CondVar_Diff_f1f2=mean(Diff_f1f2.*Diff_f1f2)-mean(Diff_f1f2)^2; %normalized by N (cf help cov)
        
    end
    orderedData(i).Data=tempData;
    
    if PLOTCONDVALUESEXTR
        % --- extrinsic -----
        % get data together 
        xdata=zeros(length(tempData),1);
        covdata=zeros(size(xdata));
        meanf1data=zeros(size(xdata));
        meanf2data=zeros(size(xdata));
        
        for run=1:length(xdata)
            xdata(run)=tempData(run).muAv;
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
        title(['Extr -> conditional variables: (Variable|mu). bin size frac=' num2str(MuBins(i))]);
        xlabel('mu')
        ylabel('mean(..|mu) resp cov(..|mu)')
        legend('Cov(field1,field2|mu)', 'E(field1|mu)', 'E(field2|mu)','Location','NW')
        grid on
    end
    if PLOTCONDVALUESINTR
        % --- intrinsic -----
        % get data together 
        xdata=zeros(length(tempData),1);
        vardata=zeros(size(xdata));
        meanDiff_f1f2data=zeros(size(xdata));
        
        for run=1:length(xdata)
            xdata(run)=tempData(run).muAv;
            vardata(run)=tempData(run).CondVar_Diff_f1f2;
            meanDiff_f1f2data(run)=tempData(run).CondMean_Diff_f1f2;
        end
        figure
        set(gcf,'WindowStyle','docked')
        clf
        plot(xdata,vardata,'.-r')
        hold on
        plot(xdata,meanDiff_f1f2data,'.-b')
        title(['Intr -> conditional variables: (Variable|mu). bin size frac=' num2str(MuBins(i))]);
        xlabel('mu')
        ylabel('mean(..|mu) resp var(..|mu)')
        legend('Var((field1-field2)|mu)', 'E((field1-field2)|mu)','Location','NW')
        grid on
        
    end
    
end % loop over all binSizes

clear Diff_f1f2 PLOTCONDVALUESEXTR PLOTCONDVALUESINTR  f1 f2 i meanDiff_f1f2data 
%clear meanf1data meanf2data n run tempData vardata xdata

%%   
% ************************************************
% CALCULATE TOTAL COVARIANCEs/EXPECTATION VALUES
% (explained_Extr+unexplained_Extr part).
% ADD VALUES TO orderedData
% ************************************************
%MuBins=[1/5,  1/10,    1/20,   1/100]; %redundant

% -----------
% law of total covariance
% ----------
% perform outer averages resp. covariances
for i=1:length(MuBins)  
    % for readability:
    tempData=orderedData(i).Data;
    % get data together
    % -- extrinsic --
    muAv_vec=[];%zeros(length(tempData),1);  % in plots named: xdata
    CondCov_field1field2_vec=[];%zeros(size(muAv));
    CondMean_field1_vec=[];%zeros(size(muAv));
    CondMean_field2_vec=[];%zeros(size(muAv));
    probMuBin_vec=[];%zeros(size(muAv));
    % -- additional intrinsic --
    CondVar_Diff_f1f2_vec=[];
    CondMean_Diff_f1f2_vec=[];
    
    % loop over each bin and extract data. ONLY USE NOT-NAN data (if data
    % contains NaN, then this bin is empty and has probability =0 -> will
    % not contribute to averages in an analytical solution)
    for n=1:length(tempData);
        if ~isnan(tempData(n).CondCov_field1field2) % bin is not empty
            %muAv_vec necessary?
            CondCov_field1field2_vec(end+1)=tempData(n).CondCov_field1field2;
            CondMean_field1_vec(end+1)=tempData(n).CondMean_field1;
            CondMean_field2_vec(end+1)=tempData(n).CondMean_field2;
            probMuBin_vec(end+1)=tempData(n).probMuBin;
            CondVar_Diff_f1f2_vec(end+1)=tempData(n).CondVar_Diff_f1f2;
            CondMean_Diff_f1f2_vec(end+1)=tempData(n).CondMean_Diff_f1f2;
        end
    end
    
    % perform averages/covariances  
    % ---- EXTRINSIC ----
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
    
    % ---- INTRINSIC ----
    %better readability
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
    
end

% some output
disp(' ')
disp('_______________________________________________________________')
disp(['Analyzed ' schnitzUseName '.'])
disp(['Inputs were extra normalized (zscore); 1=true,0=false : ' num2str(NormalizeStdDev)]);
disp(['Used fields ' field1 ' & ' field2 ])
disp(['and conditioned on ' field3 ' without history.'])
disp('  ')
disp('Extrinsic')
disp(['BinFrac     Explained_Extr  Unexplained_Extr   Total        Frac_Expl     Frac_Unexpl'])
for i=1:length(orderedData)
    disp([num2str(MuBins(i)) '         ' num2str(orderedData(i).Cov_of_Mean_Explained_Extr) '        ' ... 
        num2str(orderedData(i).Mean_of_Cov_Unexplained_Extr) ...
        '            '   num2str(orderedData(i).totalCov_Extr)  '      '  ...
        num2str(orderedData(i).Cov_of_Mean_Explained_Extr/orderedData(i).totalCov_Extr) ...
        '       ' num2str(orderedData(i).Mean_of_Cov_Unexplained_Extr/orderedData(i).totalCov_Extr)]);
end

% *******Alternative*********** BLUBB
%disp(['BinFrac     Explained_Extr  Unexplained_Extr   Total        Frac_Expl     Frac_Unexpl ALTERNATIVE'])
%for i=1:length(orderedData)
%    disp([num2str(MuBins(i)) '         ' num2str(orderedData(i).Cov_of_Mean_Expl_Extr_ALT) '        ' ... 
%        num2str(orderedData(i).Mean_of_Cov_Unexpl_Extr_ALT) ...
%        '            '   num2str(orderedData(i).totalCov_Extr)  '      '  ...
%        num2str(orderedData(i).Cov_of_Mean_Expl_Extr_ALT/orderedData(i).totalCov_Extr) ...
%        '       ' num2str(orderedData(i).Mean_of_Cov_Unexpl_Extr_ALT/orderedData(i).totalCov_Extr)]);
%end
% ******************

disp(' ')
disp('Intrinsic')
disp(['BinFrac     Explained_Intr  Unexplained_Intr   Total        Frac_Expl     Frac_Unexpl'])
for i=1:length(orderedData)
    disp([num2str(MuBins(i)) '         ' num2str(orderedData(i).Var_of_Mean_Explained_Intr) '        ' ... 
        num2str(orderedData(i).Mean_of_Var_Unexplained_Intr) ...
        '           '   num2str(orderedData(i).totalVar_Intr)  '     '  ...
        num2str(orderedData(i).Var_of_Mean_Explained_Intr/orderedData(i).totalVar_Intr) ...
        '       ' num2str(orderedData(i).Mean_of_Var_Unexplained_Intr/orderedData(i).totalVar_Intr)]);
end

% *******Alternative*********** BLUBB -> lead to exactly same result
% (implementation thus correct...)
%disp(['BinFrac     Explained_Intr  Unexplained_Intr   Total        Frac_Expl     Frac_Unexpl ATLTERNATIVE'])
%for i=1:length(orderedData)
%    disp([num2str(MuBins(i)) '         ' num2str(orderedData(i).Var_of_Mean_Expl_Intr_ALT) '        ' ... 
%        num2str(orderedData(i).Mean_of_Var_Unexpl_Intr_ALT) ...
%        '           '   num2str(orderedData(i).totalVar_Intr)  '     '  ...
%        num2str(orderedData(i).Var_of_Mean_Expl_Intr_ALT/orderedData(i).totalVar_Intr) ...
%        '       ' num2str(orderedData(i).Mean_of_Var_Unexpl_Intr_ALT/orderedData(i).totalVar_Intr)]);
%end
% ******************

disp('_______________________________________________________________')

clear CondCov_field1field2_vec CondMean_Diff_f1f2_vec CondMean_field1_vec
clear CondMean_field2_vec CondVar_Diff_f1f2_vec diff i n muAv_vec m probMuBin_vec
clear tempData v x y


%%   
% ************************************************
% INVESTIGATE DEPENDENCE OF EXPLAINED FRACTION ON #BINS
% extract data and create plot
% ************************************************

% create matrix with: fracbinSize --- #bins --- fracexplExtr --- fracexplIntr
matBinDependence=zeros(length(MuBins),4);
matBinDependence(:,1)=MuBins;
if ~NormalizeStdDev
    BinSizeAbsArray=MuBins*mumean;
    NumBinsArray=ceil(max(mymatrix(:,3))./BinSizeAbsArray);
    matBinDependence(:,2)=NumBinsArray;
end
matBinDependence(:,3)=[orderedData.Cov_of_Mean_Explained_Extr]./[orderedData.totalCov_Extr];
matBinDependence(:,4)=[(orderedData.Var_of_Mean_Explained_Intr)]./[(orderedData.totalVar_Intr)];

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
        xlim([0 40])
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
        xlim([0 20])
        xlabel('#bins')
        ylabel('frac_explained','Interpreter','None');
        set(gca, 'XDir','reverse')
        if yy(2)>0  % for 'normal' data always the case
           set(gca,'ylim',[0 yy(2)]);
        end
        grid on
    end
    annotation(gcf,'textbox',[0.191 0.485 0.82525 0.0412500000000001],...
        'String',{['Plotted ' num2str(size(mymatrix,1)) ' Schnitzes of ' schnitzUseName '.'], ...
        [field1 ', ' field2 ' conditioned on ' field3 '.']},...
        'FitBoxToText','on','Interpreter','none', 'FontSize',12);
end

clear PLOTBINDEPENDENCE matBinDependence meanf1data meanf2data run xdata vardata covdata


%%   
% ************************************************
% Extract conditional averages and variances for specific mu-bin into new
% variables
% ************************************************
CondMean_f1=[]; % conditional mean of field 1 (e.g. yfprate)
CondMean_f2=[];
CondCov_f1f2=[];
CondMean_f1minf2=[];   % conditional mean of field1 - field2   :   <yfp-cfp|mu>
CondVar_f1minf2=[];    % conditional Variance of field1 - field2   :   Var(yfp-cfp|mu)
MuVec=[];    % vector which contains the mean of each mu-bin
ProbMuBin=[];

for i=1:size(orderedData,2) % loop over mu-bin-sizes and find the correct one
    if orderedData(1,i).muBinSize==SpecialMuBinCondValues
        subData=orderedData(1,i).Data;
        for runmu=1:length(subData) % loop over each mu-bin
            CondMean_f1=[CondMean_f1, subData(runmu).CondMean_field1]; % conditional mean of field 1 (e.g. yfprate)
            CondMean_f2=[CondMean_f2, subData(runmu).CondMean_field2];
            CondCov_f1f2=[CondCov_f1f2, subData(runmu).CondCov_field1field2];
            CondMean_f1minf2=[CondMean_f1minf2, subData(runmu).CondMean_Diff_f1f2];   % conditional mean of field1 - field2   :   <yfp-cfp|mu>
            CondVar_f1minf2=[CondVar_f1minf2, subData(runmu).CondVar_Diff_f1f2];
            MuVec=[MuVec, subData(runmu).muAv];
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
    end
end

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
    title('Conditional Variance (Var(f1-f2|mu). Intr noise. no norm')
    xlabel(['mu (binfrac:' num2str(SpecialMuBinCondValues)])
    ylabel('var(...|mu)')
    grid on
end

clear i runmu idx PLOTCONDVARCOV