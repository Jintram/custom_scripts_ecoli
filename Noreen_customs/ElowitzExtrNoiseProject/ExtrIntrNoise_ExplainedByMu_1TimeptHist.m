% Calculates to which part mu-related fluctuations are responsible for
% extrinsic and intrinsic noise
% production rates (or concentrations) must be contained in field1 and
% field2 and will be normalized to average=1


% ***************
% Main Ouput: 'ordereredData' (contains all information)
% ***************

% ***************
% Nomenclautre:
% mu: growth rate at current time step (at same time as field1, field2)
% muHist: growth rate of time step before (DeltaTime=time in between fluor
%        frames)
% ***************

% ******************
% TAKES ONE TIME POINT OF HISTORY INTO ACCOUNT - altered version of
% ExtrIntrNoise_ExplainedByMu.m
% If data point is the first one in history of a specific schnitz, the last
% data point of its parent will be used
% ******************

% **********************************************
% INITIALIZE AND EXTRACT DATA
% **********************************************

% *** FACULTATIVE ***
% run initiation and selection from schnitzcells excel file
%load '\\biofysicasrv\Users2\Walker\ExperimentsCollectedData\schnitzcellsCollectionForNoiseExplByMu20120913.mat';
%

%%
% **** ADJUST *****
PLOTMUHIST=1;
PLOTCONDVALUESEXTR=1;
PLOTCONDVALUESINTR=0;

schnitzUseName='sch';  %schnitzcells selection
field1='dY5_len_smooth3_dt_subtr';  
field2='dC5_len_smooth3_dt_subtr';
field3='muP15_fitNew_subtr'; 
%field1='dY5_sum_dt_subtr';  % production rates, maybe concentrations
%field2='dC5_sum_dt_subtr';  % production rates, maybe concentrations
%field3='muP15_fitNew_subtr'; %growth rate!
NormalizeStdDev=0; % ** if =1, normalizes the standard deviations of field1 and field2 to =1. 
                   %sets mean to =0 (maybe way of normalization still needs to be improved!!)
                   % Now: ALSO field3 (MU) is standardized!!
                   % All fields get mean=0!!! (careful when binning)
                   % ** if=0, field1 and field2 are normalized to mean=1,
                   % stddev is untouched.
MuBins=[0.1 0.3 ];%1/20];  % size of one bin, in fraction of <mu>. can be array.
% *****************

% get schitzcells structure (above: only name)
eval(['schnitzUse=' schnitzUseName ';']);
% extract data
mymatrix=zeros(0,4); %  *** field1 - field2 - field3 (mu) - field3 (muHist) ***.   (can include weights here)

% loop over all schnitzes
for i=1:length(schnitzUse)
    s=schnitzUse(i);
    % use only schnitzes with useForPlot=1
    if s.useForPlot==1
        % check if fields actually contain data (length>0)
        % Note: e.g. prod.rates and mu do not always have same array length
        % (rates can be 1 shorter for late data)
        minLength=min([length(s.(field1)),length(s.(field2)),length(s.(field3))]);
        if minLength>0
            % loop over time points of schnitz s, check if values are NaN
            % and get history data point of field3 (mu). Add data to
            % mymatrix
            for time=1:minLength
                NanTestOk=isnan(s.(field1)(time)) + isnan(s.(field2)(time)) + isnan(s.(field3)(time));
                if NanTestOk==0
                    % get hist data point
                    if time==1 %first data pt of schnitz
                        if s.P~=0 % parent exists
                            parentfield3=schnitzUse(s.P).(field3); % fields should/must all have same length for parent
                            if length(parentfield3)>0 & ~isnan(parentfield3(end))
                                mymatrix=[mymatrix; s.(field1)(time), s.(field2)(time), s.(field3)(time), parentfield3(end)];
                            end
                        end
                    else
                        if ~isnan(s.(field1)(time-1))
                            mymatrix=[mymatrix; s.(field1)(time), s.(field2)(time), s.(field3)(time), s.(field3)(time-1)];
                        end
                    end
                end
            end
        end
    end
end

%normalize
mymatrix(:,1)=mymatrix(:,1)./mean(mymatrix(:,1));
mymatrix(:,2)=mymatrix(:,2)./mean(mymatrix(:,2));
if NormalizeStdDev
    mymatrix(:,1)=zscore(mymatrix(:,1));
    mymatrix(:,2)=zscore(mymatrix(:,2));
    % get mu bin sizes before subtracting mean. divide by std deviation
    % (->zscore rescaling)
    MuBinsAbsStd=(MuBins.*mean(mymatrix(:,3)))./std(mymatrix(:,3));
    mymatrix(:,3)=zscore(mymatrix(:,3));
    mymatrix(:,4)=zscore(mymatrix(:,4));
    
end
% force negative growth rates to be =0 (blubb) if mu not extra normalized. Alternatively, this data
% could also be deleted. step necessary for consistent mu-binning (don't
% forget some data)
if ~NormalizeStdDev
    mymatrix(mymatrix(:,3)<0,3)=0;
    mymatrix(mymatrix(:,4)<0,4)=0;
end
 
clear NanTestOk parentfield3 minLength s i
clear schnitzUse time

%%
% **********************************************
% CREATE STRUCT THAT CONTAINS DATA BINNED ACCORDING TO MU & MUHIST (orderedData)
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
% Similar to the non-history case, each entry tempData(k) contains all the
% information about a mu-bin. Now, each entry corresponds to a bin of mu
% a-n-d muHist. Data is stored in one long vector of structs (alternative
% would be a matrix). If n mu-bins exist, the length of the vector is n*n.
% Ordering:                 k        mu         muHist    
%          tempData(k)      1        0           0
%                           2        0           1
%                                 .....
%                           n        0          max
%                          n+1       1           0
%                                 .....
%                          n^2       max        max
% (the 0,1,2,max is not an absolute value for mu but the for simplicity just
% the bin number)
% tempData(k) contains:
% - muMin: lower border for mu for which data points land in this
%          container.
% - muMax: upper border of this bins.
% - muAv: average (centered) mu of this bin. In the middle of muMin and
%         muMax
% - muHistMin: lower border for muHist for which data points land in this
%          container.
% - muHistMax: upper border of this bins.
% - muHistAv: average (centered) muHist of this bin. In the middle of muMin and
%         muMax
% - subMatrixMuBin: submatrix of mymatrix for which mu and muHist are within the borders
% - probMuBin: fraction of all data points that lie within this muBin

% -----------------------------------------------------
% if not standardized variables (especially mu!)
% -----------------------------------------------------
NumBinsVec=zeros(size(MuBins)); %store the #of mu-bins for each fraction size
if ~NormalizeStdDev
    for i=1:length(MuBins)
        tempData=struct([]);
        % absolute size of one Bin (in dbl/h)
        BinSizeAbs=MuBins(i)*mumean;
        % get necessary number of bins to capture also highest mu. (bins start
        % at mu=0, all data with mu<0 is forced to mu=0)
        NumBins=ceil(max([mymatrix(:,3);mymatrix(:,4)])/BinSizeAbs);
        NumBinsVec(i)=NumBins;
        % loop over each mu-bin and muHist-bin and extract data that lies within this bin
        for n_mu=1:NumBins
            for n_muHist=1:NumBins
                n=(n_mu-1)*NumBins+n_muHist; % index in tempData
                tempData(n).muMin=(n_mu-1)*BinSizeAbs;
                tempData(n).muMax=n_mu*BinSizeAbs;
                tempData(n).muAv=(n_mu-0.5)*BinSizeAbs;
                tempData(n).muHistMin=(n_muHist-1)*BinSizeAbs;
                tempData(n).muHistMax=n_muHist*BinSizeAbs;
                tempData(n).muHistAv=(n_muHist-0.5)*BinSizeAbs;
                idx=find(mymatrix(:,3)>=tempData(n).muMin & mymatrix(:,3)<tempData(n).muMax ...
                     & mymatrix(:,4)>=tempData(n).muHistMin & mymatrix(:,4)<tempData(n).muHistMax);
                tempData(n).subMatrixMuBin=mymatrix(idx,:);
                tempData(n).probMuBin=length(idx)./size(mymatrix,1);
            end
        end
        orderedData(i).Data=tempData;

        if PLOTMUHIST
            % get data together
            mudata=zeros(length(tempData),1);
            muHistdata=zeros(length(tempData),1);
            probdata=zeros(length(tempData),1);
            for run=1:length(mudata)
                mudata(run)=tempData(run).muAv;
                muHistdata(run)=tempData(run).muHistAv;
                probdata(run)=tempData(run).probMuBin;
            end
            %cumulative distributions for mu and muHist (should be the
            %same)
            muBinCenters=unique(mudata);  % the same as unique(muHistdata)
            probmuCum=zeros(size(muBinCenters));
            probmuHistCum=zeros(size(muBinCenters));
            for j=1:length(muBinCenters)
                idx=find(mudata==muBinCenters(j));
                probmuCum(j)=sum(probdata(idx));
                idx=find(muHistdata==muBinCenters(j));
                probmuHistCum(j)=sum(probdata(idx));
            end
            % prob distribution in matrix form
            muprobMatrix=zeros(NumBins,NumBins);
            for mat=1:NumBins
            muprobMatrix(mat,1:NumBins)=probdata((mat-1)*NumBins+1:mat*NumBins);
            end

            figure
            set(gcf,'OuterPosition',[200 200 600 800])
            clf
            subplot(2,1,1)
            colormap ([0 0 1; 1 0 0])
            bar(muBinCenters,[probmuCum,probmuHistCum])
            legend('mu','muHist')
            title(['growth rate distribution']);
            xlabel('mu')
            ylabel('fraction')
            annotation(gcf,'textbox',...
            [0.661714285714286 0.680952380952381 0.193642857142857 0.145238095238098],...
             'String',{['<mu>=' num2str(mumean)],['bin size frac=' num2str(MuBins(i))]}, 'FitBoxToText','on');
         
            subplot(2,1,2)
            colormap('default')
            surf(muBinCenters,muBinCenters,muprobMatrix,'LineStyle','None')
            xlabel('muHist') % stimmt das?
            ylabel('mu')
            zlabel('joined prob')
            xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
            ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
            colorbar
           
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
        NumBins=ceil(max([mymatrix(:,3);mymatrix(:,4)])/BinSizeAbs-min([mymatrix(:,3);mymatrix(:,4)])/BinSizeAbs);
        NumBinsVec(i)=NumBins;
        MuStart=min([mymatrix(:,3);mymatrix(:,4)]);
        % loop over each mu-bin and muHist-bin and extract data that lies within this bin
        for n_mu=1:NumBins
            for n_muHist=1:NumBins
                n=(n_mu-1)*NumBins+n_muHist; % index in tempData
                tempData(n).muMin=MuStart+(n_mu-1)*BinSizeAbs;
                tempData(n).muMax=MuStart+n_mu*BinSizeAbs;
                tempData(n).muAv=MuStart+(n_mu-0.5)*BinSizeAbs;
                tempData(n).muHistMin=MuStart+(n_muHist-1)*BinSizeAbs;
                tempData(n).muHistMax=MuStart+n_muHist*BinSizeAbs;
                tempData(n).muHistAv=MuStart+(n_muHist-0.5)*BinSizeAbs;
                idx=find(mymatrix(:,3)>=tempData(n).muMin & mymatrix(:,3)<tempData(n).muMax ...
                     & mymatrix(:,4)>=tempData(n).muHistMin & mymatrix(:,4)<tempData(n).muHistMax);
                tempData(n).subMatrixMuBin=mymatrix(idx,:);
                tempData(n).probMuBin=length(idx)./size(mymatrix,1);
            end
        end
        orderedData(i).Data=tempData;

        if PLOTMUHIST
            % get data together
            mudata=zeros(length(tempData),1);
            muHistdata=zeros(length(tempData),1);
            probdata=zeros(length(tempData),1);
            for run=1:length(mudata)
                mudata(run)=tempData(run).muAv;
                muHistdata(run)=tempData(run).muHistAv;
                probdata(run)=tempData(run).probMuBin;
            end
            %cumulative distributions for mu and muHist (should be the
            %same)
            muBinCenters=unique(mudata);  % the same as unique(muHistdata)
            probmuCum=zeros(size(muBinCenters));
            probmuHistCum=zeros(size(muBinCenters));
            for j=1:length(muBinCenters)
                idx=find(mudata==muBinCenters(j));
                probmuCum(j)=sum(probdata(idx));
                idx=find(muHistdata==muBinCenters(j));
                probmuHistCum(j)=sum(probdata(idx));
            end
            % prob distribution in matrix form
            muprobMatrix=zeros(NumBins,NumBins);
            for mat=1:NumBins
            muprobMatrix(mat,1:NumBins)=probdata((mat-1)*NumBins+1:mat*NumBins);
            end

            figure
            set(gcf,'OuterPosition',[200 200 600 800])
            clf
            subplot(2,1,1)
            colormap ([0 0 1; 1 0 0])
            bar(muBinCenters,[probmuCum,probmuHistCum])
            legend('mu','muHist')
            title(['growth rate distribution']);
            xlabel('mu')
            ylabel('fraction')
            annotation(gcf,'textbox',...
            [0.661714285714286 0.680952380952381 0.193642857142857 0.145238095238098],...
             'String',{['<mu>=' num2str(mumean)],['bin size frac=' num2str(MuBins(i))]}, 'FitBoxToText','on');
         
            subplot(2,1,2)
            colormap('default')
            surf(muBinCenters,muBinCenters,muprobMatrix,'LineStyle','None')
            xlabel('muHist') % stimmt das?
            ylabel('mu')
            zlabel('joined prob')
            xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
            ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
            colorbar
            
        end
    end % loop over all binSizes
    
end     

clear BinSizeAbs PLOTMUHIST i idx n run tempData NumBins j mat 
clear muBinCenters muHistdata mudata muprobMatrix n_mu n_muHist probdata 
clear probmuCum probmuHistCum 
        
        
%%
% **********************************************
% CALCULATE CONDITIONAL (CO)VARIANCES AND EXPECTATION VALUES
% **********************************************   
% i.e.: calculate them seperately for each mu/muHist-Bin (and of course seperately
% for different binning sizes)
%PLOTCONDVALUESEXTR=0;
%PLOTCONDVALUESINTR=0;


% -----------
% conditional covariance
% ----------
%loop over all binning sizes (i)
for i=1:length(MuBins)  
    NumBins=NumBinsVec(i);
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
        %  *** get data together ***
        mudata=zeros(length(tempData),1);
        muHistdata=zeros(size(mudata));
        covdata=zeros(size(mudata));
        meanf1data=zeros(size(mudata));
        meanf2data=zeros(size(mudata));
        probdata=zeros(size(mudata));
        for run=1:length(mudata)
            mudata(run)=tempData(run).muAv;
            muHistdata(run)=tempData(run).muHistAv;
            covdata(run)=tempData(run).CondCov_field1field2;
            meanf1data(run)=tempData(run).CondMean_field1;
            meanf2data(run)=tempData(run).CondMean_field2;
            probdata(run)=tempData(run).probMuBin;
        end
        muBinCenters=unique(mudata);  % the same as unique(muHistdata)
        % prob distribution and cov/mean in matrix form
        muprobMatrix=zeros(NumBins,NumBins);
        meanf1Matrix=zeros(size(muprobMatrix));
        meanf2Matrix=zeros(size(muprobMatrix));
        covMatrix=zeros(size(muprobMatrix));
        for mat=1:NumBins
            % a row corresponds to a specific mu, a column to a specific
            % muHist
            muprobMatrix(mat,1:NumBins)=probdata((mat-1)*NumBins+1:mat*NumBins);
            meanf1Matrix(mat,1:NumBins)=meanf1data((mat-1)*NumBins+1:mat*NumBins);
            meanf2Matrix(mat,1:NumBins)=meanf2data((mat-1)*NumBins+1:mat*NumBins);
            covMatrix(mat,1:NumBins)=covdata((mat-1)*NumBins+1:mat*NumBins);
        end
        %get color ranges for plots ***
        colorcutoff=3; % =1 colorrange from highest to lowest element. =2 range from 2nd highest to 2nd lowest etc
        % ***
        idx=~isnan(meanf1data);
        meanf1datasorted=sort(meanf1data(idx));
        meanf1color=[meanf1datasorted(colorcutoff), meanf1datasorted(end+1-colorcutoff)];
        idx=~isnan(meanf2data);
        meanf2datasorted=sort(meanf2data(idx));
        meanf2color=[meanf2datasorted(colorcutoff), meanf2datasorted(end+1-colorcutoff)];
        idx=~isnan(covdata);
        covdatasorted=sort(covdata(idx));
        covcolor=[covdatasorted(colorcutoff), covdatasorted(end+1-colorcutoff)];
        
        % *****************************************
        % *** plot figure extr ***
        figure
        clf
        hold on
        set(gcf,'OuterPosition',[50 50 1600 800]);
        
        subplot(2,4,8,'Position',[0.78 0.05 0.19 0.38])   % probs. images
        h=imagesc(muBinCenters,muBinCenters, muprobMatrix);
        set(h,'alphadata',~isnan(muprobMatrix));
        set(gca,'YDir','normal');
        xlabel('muHist')
        ylabel('mu')
        %xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        %ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])      
        colormap hot
        freezeColors % freezes colormap
        hc=cbfreeze; % creates colorbar and freezes colors
        ylabel(hc,'prob to be in certain bin','FontSize',12)
        
        colormap jet
        
        subplot(2,4,1,'Position',[0.03 0.55 0.19 0.38])  % 1st plot. surf   
        surf(muBinCenters,muBinCenters, muprobMatrix, covMatrix,'LineStyle','None')
        view([0 90]) % set direction of view -> in x/y plane
        xlabel('muHist')
        ylabel('mu')
        xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'cov(f1,f2|mu)','FontSize',12)
        caxis(covcolor)
        
        subplot(2,4,5,'Position',[0.03 0.05 0.19 0.38])   % 1st plot. imagesc  
        h=imagesc(muBinCenters,muBinCenters, covMatrix);
        set(h,'alphadata',~isnan(covMatrix));
        set(gca,'YDir','normal');
        xlabel('muHist')
        ylabel('mu')
        %xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        %ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'cov(f1,f2|mu)','FontSize',12)
        caxis(covcolor)
        
        
        subplot(2,4,2,'Position',[0.28 0.55 0.19 0.38])   % 2nd plot. surf  
        surf(muBinCenters,muBinCenters, muprobMatrix, meanf1Matrix,'LineStyle','None')
        view([0 90]) % set direction of view -> in x/y plane
        xlabel('muHist')
        ylabel('mu')
        xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'mean(f1|mu)','FontSize',12)
        caxis(meanf1color)
        %general title (centered)
        title(['Extr -> conditional variables: (Variable|mu). bin size frac=' num2str(MuBins(i))],'FontSize',14);
        
        subplot(2,4,6,'Position',[0.28 0.05 0.19 0.38])   % 2nd plot. imagesc  
        h=imagesc(muBinCenters,muBinCenters, meanf1Matrix);
        set(h,'alphadata',~isnan(meanf1Matrix));
        set(gca,'YDir','normal');
        xlabel('muHist')
        ylabel('mu')
        %xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        %ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'mean(f1|mu)','FontSize',12)
        caxis(meanf1color)
        title('Correct bin borders and color below.','FontSize',14);
        
        subplot(2,4,3,'Position',[0.53 0.55 0.19 0.38])   % 3rd plot. surf       
        surf(muBinCenters,muBinCenters, muprobMatrix, meanf2Matrix,'LineStyle','None')
        view([0 90]) % set direction of view -> in x/y plane
        xlabel('muHist')
        ylabel('mu')
        xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'mean(f2|mu)','FontSize',12)
        caxis(meanf2color)
        
        subplot(2,4,7,'Position',[0.53 0.05 0.19 0.38])   % 3rd plot. imagesc  
        h=imagesc(muBinCenters,muBinCenters, meanf2Matrix);
        set(h,'alphadata',~isnan(meanf2Matrix));
        set(gca,'YDir','normal');
        xlabel('muHist')
        ylabel('mu')
        %xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        %ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'mean(f2|mu)','FontSize',12)
        caxis(meanf2color)
        
        annotation(gcf,'textbox',...
        [0.802767676767677 0.728813559322034 0.110479797979798 0.0395480225988701],...
        'String',{'upper row: 3dim surf plots'});
        % *****************************************
    end
    
    if PLOTCONDVALUESINTR
        % --- intrinsic -----
        % *** get data together ***
        mudata=zeros(length(tempData),1);
        muHistdata=zeros(size(mudata));
        vardata=zeros(size(mudata));
        meanDiff_f1f2data=zeros(size(mudata));
        probdata=zeros(size(mudata));
        for run=1:length(mudata)
            mudata(run)=tempData(run).muAv;
            muHistdata(run)=tempData(run).muHistAv;
            vardata(run)=tempData(run).CondVar_Diff_f1f2;
            meanDiff_f1f2data(run)=tempData(run).CondMean_Diff_f1f2;
            probdata(run)=tempData(run).probMuBin;
        end
        
        % prob distribution and cov/mean in matrix form
        muprobMatrix=zeros(NumBins,NumBins);
        varMatrix=zeros(size(muprobMatrix));
        meanDiff_f1f2dataMatrix=zeros(size(muprobMatrix));
        for mat=1:NumBins
            % a row corresponds to a specific mu, a column to a specific
            % muHist
            muprobMatrix(mat,1:NumBins)=probdata((mat-1)*NumBins+1:mat*NumBins);
            varMatrix(mat,1:NumBins)=vardata((mat-1)*NumBins+1:mat*NumBins);
            meanDiff_f1f2dataMatrix(mat,1:NumBins)=meanDiff_f1f2data((mat-1)*NumBins+1:mat*NumBins);
        end
        %get color ranges for plots ***
        colorcutoff=3; % =1 colorrange from highest to lowest element. =2 range from 2nd highest to 2nd lowest etc
        % ***
        idx=~isnan(vardata);
        vardatasorted=sort(vardata(idx));
        varcolor=[vardatasorted(colorcutoff), vardatasorted(end+1-colorcutoff)];
        idx=~isnan(meanDiff_f1f2data);
        meanDiff_f1f2datasorted=sort(meanDiff_f1f2data(idx));
        meanDiff_f1f2color=[meanDiff_f1f2datasorted(colorcutoff), meanDiff_f1f2datasorted(end+1-colorcutoff)];
        
        % *****************************************
        % *** plot figure intr ***
        figure
        clf
        hold on
        set(gcf,'OuterPosition',[50 50 1400 1000]);
        
        subplot(2,3,6)   % probs. imagesc
        h=imagesc(muBinCenters,muBinCenters, muprobMatrix);
        set(h,'alphadata',~isnan(muprobMatrix));
        set(gca,'YDir','normal');
        xlabel('muHist')
        ylabel('mu')
        %xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        %ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])      
        colormap hot
        freezeColors % freezes colormap
        hc=cbfreeze; % creates colorbar and freezes colors
        ylabel(hc,'prob to be in certain bin','FontSize',12)
        
        colormap jet
        
        
        subplot(2,3,1)   % 1st plot. surf
        surf(muBinCenters,muBinCenters, muprobMatrix, varMatrix,'LineStyle','None')
        view([0 90]) % set direction of view -> in x/y plane
        xlabel('muHist')
        ylabel('mu')
        xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'var(f1-f2|mu)','FontSize',12)
        caxis(varcolor)
       
        subplot(2,3,4)   % 1st plot. imagesc  
        h=imagesc(muBinCenters,muBinCenters, varMatrix);
        set(h,'alphadata',~isnan(varMatrix));
        set(gca,'YDir','normal');
        xlabel('muHist')
        ylabel('mu')
        %xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        %ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'var(f1-f2|mu)','FontSize',12)
        caxis(varcolor)
        
        subplot(2,3,2)   % 2nd plot. surf   
        surf(muBinCenters,muBinCenters, muprobMatrix, meanDiff_f1f2dataMatrix,'LineStyle','None')
        view([0 90]) % set direction of view -> in x/y plane
        xlabel('muHist')
        ylabel('mu')
        xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'mean(f1-f2|mu)','FontSize',12)
        caxis(meanDiff_f1f2color)
        %general title
        title(['Intr -> conditional variables: (Variable|mu). bin size frac=' num2str(MuBins(i))],'FontSize',14);
 
        subplot(2,3,5)   % 2nd plot. imagesc  
        h=imagesc(muBinCenters,muBinCenters, meanDiff_f1f2dataMatrix);
        set(h,'alphadata',~isnan(meanDiff_f1f2dataMatrix));
        set(gca,'YDir','normal');
        xlabel('muHist')
        ylabel('mu')
        %xlim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        %ylim([1.5*muBinCenters(1)-0.5*muBinCenters(2) 1.5*muBinCenters(end)-0.5*(muBinCenters(end-1))])
        colorbar
        ylabel(colorbar,'mean(f1-f2|mu)','FontSize',12)
        caxis(meanDiff_f1f2color)
        %general title
        title('Correct bin borders and color below.','FontSize',14);
        
        annotation(gcf,'textbox',...
        [0.702767676767677 0.728813559322034 0.110479797979798 0.0395480225988701],...
        'String',{'upper row: 3dim surf plots'});
        % *****************************************
    end
    
end % loop over all binSizes

clear Diff_f1f2 PLOTCONDVALUESEXTR PLOTCONDVALUESINTR covdata f1 f2 i meanDiff_f1f2data 
clear meanf1data meanf2data n run tempData vardata mudata NumBins colorcutoff covMatrix
clear covdatasorted h hc idx mat meanDiff_f1f2color meanDiff_f1f2dataMatrix
clear meanDiff_f1f2datasorted meanf1Matrix meanf1color meanf1datasorted
clear meanf2Matrix meanf2color meanf2datasorted muBinCenters muHistdata muprobMatrix
clear probdata varMatrix varcolor vardatasorted covcolor

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
    muAv_vec=[];%zeros(length(tempData),1);  % in plots named: mudata
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
    
    % perform averages/covariances  %BLUBB NOCHMAL GENAU KONTROLLIEREN!!!
    % HOW TO CHECK/BACKUP?
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
    
end
    
% some output
disp(' ')
disp('_______________________________________________________________')
disp(['Analyzed ' schnitzUseName '.'])
disp(['Inputs were extra normalized (zscore); 1=true,0=false : ' num2str(NormalizeStdDev)]);
disp(['Used fields ' field1 ' & ' field2 ])
disp(['and conditioned on ' field3 ' with 1 history time point (2 points in total).'])
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
disp('_______________________________________________________________')

clear CondCov_field1field2_vec CondMean_Diff_f1f2_vec CondMean_field1_vec
clear CondMean_field2_vec CondVar_Diff_f1f2_vec diff i n muAv_vec m probMuBin_vec
clear tempData v x y