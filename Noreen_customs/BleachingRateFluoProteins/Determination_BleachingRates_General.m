% --------------------------------------------------------------------
% DETERMINATION OF BLEACHING RATE FOR AN ARBITRARY COLOR, BASED ON 
% A SCHNITZCELLS STRUCTURE.
% This script is used in the context of fluorescent protein characterization.
% --------------------------------------------------------------------

% ***************************************************************
% IMPORTANT: FOR CORRECT UNEVEN-ILLUMINATION-CORRECTION (=need location of
% cell in full image) IT IS ESSENTIAL THAT IMAGE IS
% - NOT CROPPED
% - 'USEFULLIMAGE' SET
% ***************************************************************


% Origin: 
% The script is strongly based on "Determination_BleachingRates_GFPmCherry.m", 
%   however * generalized to any color
%           * specific to one color per experimental dataset (double labels
%             need to be run twice)
%           * takes uneven illumination into account (cells at edge are
%             bleached less. Larger new camera chip!)
% 
% Creation Date: NW 2015-07-07
% --------------------------------------------------------------------

% ------------
% EXPERIMENT:
% ------------
% Cells are treated with antibiotics, protein production stops. A
% time-lapse movie is acquired, fluo proteins illuminated (at every phase 
% image) and thus bleached.
% Example datasets: 2015-06-11, 2015-06-09Bleach

% ------------
% INPUT DATA NEEDED
% ------------
% * schnitzcells data structure, full tracking (growth, production rate, ... not needed). 
%   name of structure: "schnitzcells_rm", with set "useForPlot".
%   Example minimal analysis template: 
%   /2015-06-11/schnitzcell_bleachingAnalysis_Cyan_2015-06-11_pos2.xlsx
%   If merged schnitzcells file: *** all schnitzes must have same # frames!!!***
% * illumination time
% * color
% * shading image

% ------------
% ANALYSIS PERFORMED
% ------------
% * For each schnitz (=each lineage since non-dividing) of a schnitzcells 
%   dataset the bleaching rate [per time] is determined. Base=2
% * Bleaching is fitted with an ***exponential*** TO BE TESTED WHETHER
%   BLEACHING IS REALLY EXPONENTIAL (CAN BE EITHER WAY)
% * Comparison of yes/no correction of uneven illumination (less light at
%   corners...): A posterioi conclusion: Correction is absolutely
%   essential.
%   Renormalized such that bleaching/illum-sec is correct in center
%   quadrant
% * Comparison of G6mean,G5mean,G5sum.

% ------------
% OUTPUT
% ------------
% ** all of the following data is corrected for uneven illumination (larger
%    bleaching in center of chip): They report on bleaching/sec in the center quadrant **
% Bleachrates_effec:        vector with bleach rates (base=2) for each schnitz.
% T12_effec:                fluo half-time (1/Bleachrates_effec)
% Bleachrates_effec_mean:   mean bleach rate
% Bleachrates_effec_std:    stddev of bleach rate
% T12_effec_mean:           mean half-time
% T12_effec_std:            stddev of half-time
% Offsets_effec:            vector with y-axis intercepts of exp fits
% Offsets_effec_mean:       mean y-axis intercept
% Offsets_effec_std:        stddev of y-axis intercept
% ** all above quantities exist also without the "effec" and are then not
% illumination corrected -> these values are only for testing **

% MeanFitBleachCurve        Fitted bleaching curve over selected time range


% ********************************************
% STEP BY STEP: WHAT TO DO
% ********************************************
% 1) requires manual adjustement (dataset etc).
% 2) optional adjustment: ad hoc remove schnitzcells from analysis (useForPlot)
% 4) optional adjustment (X6mean vs X5sum etc, plotting options)
% 

% DIFFERENT CELLS:
% 1) load schnitzcells data, shading image and specify illumination times,
%     and color
% 2) extract fluorescence-over-time matrices from schnitzcell vectors.
%    normalize intensity.
%    Correct for uneven illumination
% 3) plot traces
% 4) determine exponential(!) bleaching rate per second. plot bleaching incl fits. determine mean and std.

%% -----------------------------------------------------------------------
% (1) LOAD SCHNITZCELLS. SPECIFY PATH AND ILLUMINATION TIME: 
% ------------------------------------------------------------------------
% ************ADJUST****************
% path (struct file in .mat must be named schnitzcells_rm)
myfile='\\biofysicasrv\Users2\Walker\TechnicalAnalysis_Calibration\BleachingRates_Main\Schnitzfiles\schnitzcells_2015-06-08pos3andpos6_bleach_mVenus.mat';
% shading image (color spec!)
myshadingfile='D:\SchnitzcellsVersions\Noreen_Develop\fluo_correction_images\Correction_microscope1-yfp-20150527-20ms.mat';
%fluo settings
illumtime=0.200; %[sec!!!]
mycolor='y';  %color (for file names)
fluoname='mVenus';
% image range to use (determine after first plot)
useimages=[1:251]; % counter starts at =1 (also if e.g pos1-..-006 is first image)
                    % 2015-06-11pos5 (mKate2): [20:350] (initial increase:
                    %                       focus or photoconversion?
% smooth traces with default span=5 (default=1)
SMOOTHTRACES=1;
% **********************************

% **** Usually don't adjust ***
load(myfile);
load(myshadingfile)
clear flatfield replace
% ok, if name was wrong: adjust schnitz file name
myschnitzcells=schnitzcells_rm;



%% --------------------------------------------------------------------------
% (2) EXTRACT FLUORESCENCE-OVER-TIME MATRICES FROM SCHNITZCELLS FILES,
%   USEFORPLOT AND EFFECTIVE ILLUMINATION STRENGTH. NORMALIZE FLUORESCENCE
% ---------------------------------------------------------------------------

numschnitzes=length(myschnitzcells);

% ------------------------------------------------------------
% Fluorescence over time
% ------------------------------------------------------------
% structure, e.g. matrix X6_mean for gfp label
%            schnitz1  -  schnitz2  -  schnitz3  -  schnitz4  -  ...
% 1st row:     G6_t0        G6_t0        G6_t0        G6_t0
% 2nd row:     G6_t1        G6_t1        G6_t1        G6_t1
% 3rd row:      ...         ....          ...          ...

% get general variable names
X6mean_name=[upper(mycolor) '6_mean'];
X5mean_name=[upper(mycolor) '5_mean'];
X5sum_name=[upper(mycolor) '5_sum'];

% allocate empty matrices
mysize=[length(myschnitzcells(1).(X6mean_name)),length(myschnitzcells)];
X6mean=zeros(mysize);   % e.g. R6_mean (conc)
X5mean=zeros(mysize);    % e.g. R5_mean (conc)
X5sum=zeros(mysize); % e.g. R5_sum (total fluo)

% fill matrices with data
for i=1:length(myschnitzcells)
    X6mean(:,i)=myschnitzcells(i).(X6mean_name)';
    X5mean(:,i)=myschnitzcells(i).(X5mean_name)';
    X5sum(:,i)=myschnitzcells(i).(X5sum_name)';
    if SMOOTHTRACES  % default span 5
        X6mean(:,i)=smooth(X6mean(:,i));
        X5mean(:,i)=smooth(X5mean(:,i));
        X5sum(:,i)=smooth(X5sum(:,i));
    end
     
        
end


% ------------------------------------------------------------
% useForPlot
% ------------------------------------------------------------
useforplot_vec=[myschnitzcells.useForPlot];
% potentially here: remove schnitzcells with bad traces adhoc -> set =0
% useforplot_vec(ii)=0; %blubb

% ------------------------------------------------------------
% Effective illumination strength
% ------------------------------------------------------------
% rescale factor for relative illumination strength depending on position.
% center cells are illumination stronger than corner cells, for same illum
% time.
% normalization is chosen such that for center quadrant cells the renormalization is ~1

% prefactors for effective illumination time: =1. "real"& "effective" the
% same. >1: effectively illuminated more (very centered cells). <1 effective illum
% is shorter (corner cells). To be multiplied with illum time
prefac_illumtime=ones(length(myschnitzcells),1);

% normalization of intensity: choose such that center quadrant has mean=1
% (actual illumination time and rescaled one are the same). why: standard
% colony is usually within this quadrant. we want to know the bleaching
% rate for these cells, with the measured (real) illumtime.

%get mean value of center quadrant of shading
numrow_shad=size(shading,1);
numcol_shad=size(shading,2);
centerquad_row=[round(numrow_shad/4) round(3*numrow_shad/4)];
centerquad_col=[round(numcol_shad/4) round(3*numcol_shad/4)];
meanshad_centerquad=mean(mean(shading(centerquad_row(1):centerquad_row(2),centerquad_col(1):centerquad_col(2))));

% normalize shading
shadingnorm=shading/meanshad_centerquad; 

% loop over all schnitzes and get relative illumination intensity
clear pos_row pos_col i dummyimage shading_schnitz shading_schnitz_mean
for i=1:length(myschnitzcells)
    pos_row=round(mean(myschnitzcells(i).ceny)); % row is the y-position! average y-pos of all frames
    pos_col=round(mean(myschnitzcells(i).cenx));
    % create a circle of radius=10 at the cell position and average
    % shading within that area
    dummyimage=zeros(size(shading));
    dummyimage(pos_row,pos_col)=1;
    dummyimage=imdilate(dummyimage, strel('disk',10));
    % shading_schnitz=dummyimage.*shadingnorm; %not needed but nice for debugging -> imagesc
    shading_schnitz_mean=mean2(shadingnorm(dummyimage>0));
    prefac_illumtime(i)=shading_schnitz_mean;
end

% -----------------
% vector with cumulative illumtime (not intensity corrected)
% -----------------
illumtime_vec=zeros(length(myschnitzcells(1).(X6mean_name)),1);
for i=1:length(illumtime_vec);
    illumtime_vec(i)=i*illumtime;
end

% ---------------------
% matrix with effective cumulative illum time
% ---------------------
%            schnitz1  -  schnitz2  -  schnitz3  -  schnitz4  -  ...
% 1st row:     1 image     1 image       ...
% 2nd row:     2 images    2 images      ...
% 3rd row:      ...         ....          ...          ...
illumtime_mat_effec=zeros(length(illumtime_vec),length(myschnitzcells));
% loop over all frames (#illum)
for i=1:length(illumtime_vec)
    for schn=1:numschnitzes
        illumtime_mat_effec(i,schn)=i*illumtime*prefac_illumtime(schn);
    end
end


% ----------------------------------------------------------
% NORMALIZE FLUORESCENCE. RESTRICT TO USED IMAGE RANGE
% ----------------------------------------------------------
%normalize fluo to value at first image=1
X6mean_norm=zeros(size(X6mean));
X5mean_norm=zeros(size(X5mean));
X5sum_norm=zeros(size(X5sum));
for i=1:numschnitzes
    X6mean_norm(:,i)=X6mean(:,i)/X6mean(1,i);
    X5mean_norm(:,i)=X5mean(:,i)/X5mean(1,i);
    X5sum_norm(:,i)=X5sum(:,i)/X5sum(1,i);
end

% restrict to used image range
illumtime_vec=illumtime_vec(useimages);
illumtime_mat_effec=illumtime_mat_effec(useimages,:);

X6mean_norm=X6mean_norm(useimages,:);
X5mean_norm=X5mean_norm(useimages,:);
X5sum_norm=X5sum_norm(useimages,:);

%% --------------------------------------------------------------
% (3) FIGURES OF TRACES
% ---------------------------------------------------------------
RANDOMCOLOR=0;
if RANDOMCOLOR  % choose random colors
    colormatrix=rand(numschnitzes,3);
else
    colormatrix=prefac_illumtime*[1 1 1];
    colormatrix=0.8*colormatrix/max2(colormatrix);
    %colormatrix(:,2)=0.5;
end

figure(1) % not illum corrected
clf
hold on
for i=1:numschnitzes
    if useforplot_vec(i)==1
        plot(illumtime_vec,X5sum_norm(:,i),'Color',colormatrix(i,:))
    end
end
% plot(X5sum_norm) % multicolor
xlabel('cum. illumtime [sec]')
ylabel([fluoname ': ' X5sum_name],'Interpreter','None')
title('not normalized for illum strength')

figure(2) % yes illum corrected
clf
hold on
for i=1:numschnitzes
    if useforplot_vec(i)==1
        plot(illumtime_mat_effec(:,i),X5sum_norm(:,i),'Color',colormatrix(i,:))
    end
end
% plot(X5sum_norm) % multicolor
xlabel('cum. illumtime [sec]')
ylabel([fluoname ': ' X5sum_name],'Interpreter','None')
title('IS normalized for illum strength')

figure(10) % not illum corrected, not intensity normalized
clf
hold on
for i=1:numschnitzes
     if useforplot_vec(i)==1
        plot(illumtime_vec,X5sum(useimages,i),'Color',colormatrix(i,:))
     end
end
% plot(X5sum_norm) % multicolor
xlabel('cum. illumtime [sec]')
ylabel([fluoname ': ' X5sum_name],'Interpreter','None')
title('NOT normalized for init fluo value & illum strength')


%% ----------------------------------------------------------------------------
% (4) DETERMINE BLEACH RATE PER SECOND. PLOT.
% ----------------------------------------------------------------------------
% Analysis performed for 1) original fluo time traces
%                        2) rescaled with effective illum time traces
% Analysis is performed for either of X6mean,X5mean,X5sum. ***ADJUST***
%       [BLUBB all fluo values produce the same bleach rates....]

% This script is based on 'CompleteScriptMaturationTimes', subsection 'Bleaching'

% *************** ADJUST ****************
% plot fitted bleach curve for each cell.
PLOTALLTRACES=0;
% fluo quantity
Fluo_norm=X5sum_norm;   % X5sum_norm, X5mean_norm, X6mean_norm
Fluo_name=X5sum_name;  % ADJUST TO LINE ABOVE!!!!
% ***************************************
Fluo_norm_effec=Fluo_norm;

mytitle=['Protein: ' fluoname  '(' mycolor ')'];

    
% % *** define a long list of colors for plotting traces ***
%mycoloror=[1 1 0.9; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; ... 
%    0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; ...
%    1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; ...
%    0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];

% vector to store bleach rates
Bleachrates=NaN(numschnitzes,1); % no account for illum intensity.
                                       % removed schnitzes: 'NaN' stays
Bleachrates_effec=NaN(numschnitzes,1); % yes account for illum intensity

Offsets=NaN(numschnitzes,1); % y-axis offset. no account for illum intensity.
Offsets_effec=NaN(numschnitzes,1); % y-axis offset. yes account for illum intensity

% ----------------------------------------
% Bleachrates1: Not normalized for illum intensity [only for reference -
% don't actually use]
% ----------------------------------------
% loop over all schnitzes
for i=1:numschnitzes
    if useforplot_vec(i)==1
        curr_fluovec=Fluo_norm(:,i); % fluo for all frames
        curr_illumvec=illumtime_vec; % cum. illumtime for all frames

        % fit an exponential bleaching function (adapted from curve fitting session & maturation time script)
        ok2 = isfinite(curr_illumvec) & isfinite(curr_fluovec);
        if ~all( ok2 )
            warning( 'Bleach factor calculation: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
        end
        % suitable start paramters for [a b]
        st2= [ 400000 -0.01];
        ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
        % Fit this model using new data
        cfbleach = fit(curr_illumvec(ok2),curr_fluovec(ok2),ft2,'Startpoint',st2);

        % get fitted coefficients
        a=cfbleach.a; bleachrate_natlog=-cfbleach.b;
        % change log_e to log_2
        bleachrate=bleachrate_natlog/log(2); % value is larger

        %optional traces poltting
      %  if PLOTALLTRACES
      %      yaxisfitbleach=a*exp(-bleachrate*curr_illumvec);
      %      figure
      %      clf 
      %      hold on
      %      set(gcf,'WindowStyle','docked')
      %      h1=plot(curr_illumvec,curr_fluovec,'.','MarkerSize',15);
      %      h2=plot(curr_illumvec,yaxisfitbleach,'r','LineWidth',2);
      %      xlabel('cum. illumtime [sec]')
      %      ylabel('fluorescence')
      %      legend(h2,['b=-' num2str(bleachrate)])
      %      title(['Bleach curve: ' mytitle  ' schnitz nr ' num2str(i) '. not norm. for intensity.'],'Interpreter','None');
      %  end

        % store all bleachrates into one vector
        Bleachrates(i)=bleachrate;
        Offsets(i)=a;
    end
end


% ----------------------------------------
% Bleachrates2: Yes normalized for illum intensity (almost same script as
% above)
% ----------------------------------------
% loop over all schnitzes
for i=1:numschnitzes
    if useforplot_vec(i)==1
        curr_fluovec=Fluo_norm(:,i); % fluo for all frames
        curr_illumvec=illumtime_mat_effec(:,i); % cum. illumtime for all frames

        % fit an exponential bleaching function (adapted from curve fitting session & maturation time script)
        ok2 = isfinite(curr_illumvec) & isfinite(curr_fluovec);
        if ~all( ok2 )
            warning( 'Bleach factor calculation: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
        end
        % suitable start paramters for [a b]
        st2= [ 400000 -0.01];
        ft2 = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
        % Fit this model using new data
        cfbleach = fit(curr_illumvec(ok2),curr_fluovec(ok2),ft2,'Startpoint',st2);

        % get fitted coefficients
        a=cfbleach.a; bleachrate_natlog=-cfbleach.b;
        % change log_e to log_2
        bleachrate=bleachrate_natlog/log(2); % value is larger

        %optional traces poltting
        if PLOTALLTRACES
            yaxisfitbleach=a*2.^(-bleachrate*curr_illumvec);
            figure
            clf 
            hold on
            set(gcf,'WindowStyle','docked')
            h1=plot(curr_illumvec,curr_fluovec,'.','MarkerSize',15);
            h2=plot(curr_illumvec,yaxisfitbleach,'r','LineWidth',2);
            xlabel('cum. illumtime [sec]')
            ylabel('fluorescence')
            legend(h2,['b=-' num2str(bleachrate)])
            title(['Bleach curve: ' mytitle  ' schnitz nr ' num2str(i) '. IS norm. for intensity.'],'Interpreter','None');
        end

        % store all bleachrates into one vector
        Bleachrates_effec(i)=bleachrate;
        Offsets_effec(i)=a;
    end
end

% ------------------------------
% Get Mean & Stdev of BleachRates. Get Half-time [sec]
% ------------------------------
idx=find(~isnan(Bleachrates));
Bleachrate_mean=mean(Bleachrates(idx)); % no illum correction
Bleachrate_std=std(Bleachrates(idx));
Bleachrate_effec_mean=mean(Bleachrates_effec(idx));
Bleachrate_effec_std=std(Bleachrates_effec(idx));

T12=1./Bleachrates;
T12_effec=1./Bleachrates_effec;

T12_mean=mean(T12(idx)); % no illum correction
T12_std=std(T12(idx));
T12_effec_mean=mean(T12_effec(idx));
T12_effec_std=std(T12_effec(idx));

Offsets_effec_mean=mean(Offsets_effec(idx));
Offsets_effec_std=std(Offsets_effec(idx));
Offsets_mean=mean(Offsets(idx));
Offsets_std=std(Offsets(idx));



% figure compare different Fluo quantitites -> bleach rate. outdated
% figure(5)
%clf
%hold on
%bar(avGFPbleach, 'FaceColor', [0.3 1 0.3])  % Lighter so error bars show up
%errorbar(avGFPbleach,stdGFPbleach, 'ks','LineWidth',2);            % Error bars use black squares
%set(gca, 'XTick', 1:3, 'XTickLabel', {'G6mean' 'G5mean'  'G5total'}) % Set ticks and tick labels
%ylabel('bleachrate (mean,std)')
%title([mytitle1 'GFP. exp bleach rate per sec'],'Interpreter','None')
%legend('Mean (SD error bars)', 'Location', 'Northwest') % Put in lower right
%box on                                         % Force box around axes
%hold off


% -----------------------------------
% FINAL PLOT: ALL TRACES & FITTED CURVE
% -----------------------------------
% Partly copied form (3)

% Fitted average bleaching curve
IllumTimeMeanCurve=illumtime_vec; % the real illumination time. all values are normalized to this time
MeanFitBleachCurve=Offsets_effec_mean*2.^(-Bleachrate_effec_mean*IllumTimeMeanCurve);

RANDOMCOLOR=0;
if RANDOMCOLOR  % choose random colors
    colormatrix=rand(numschnitzes,3);
else
    colormatrix=prefac_illumtime*[1 1 1];
    colormatrix=0.8*colormatrix/max2(colormatrix);
    %colormatrix(:,2)=0.5;
end

figure(3) % yes illum corrected
clf
hold on
for i=1:numschnitzes
    if useforplot_vec(i)==1
        plot(illumtime_mat_effec(:,i),X5sum_norm(:,i),'Color',colormatrix(i,:))
    end
end
plot(IllumTimeMeanCurve,MeanFitBleachCurve,'r','LineWidth',2)
% plot(X5sum_norm) % multicolor
xlabel('cum. illumtime [sec]')
ylabel([fluoname ': ' X5sum_name],'Interpreter','None')
legend(['T12=' num2str(T12_effec_mean)  '+-' num2str(T12_effec_std) 'sec'] )
title('IS normalized for illum strength')



% ------------------------------------
% SOME OUTPUT
% ------------------------------------
disp('  ')
disp('--------------------------------------------------------')
disp('DETERMINATION OF BLEACHING RATE')
disp('uneven illumination corrected. bleachrate = bleaching/illum-sec in center quadrant)')
disp('Fitted function: fluovalue=offset*2^(-bleachrate*total_illum)')
disp('--------------------------------------------------------')
disp(['Analyzed dataset ' myfile ' .'])
disp(['Fluorophore: ' fluoname '(' mycolor ')'])
disp(['Fluoparameter: ' Fluo_name ', illumtime per frame: ' num2str(illumtime) ... 
    'sec, #frames used: ' num2str(length(useimages)) '.'])
disp(['# schnitzes (traces) used: ' num2str(sum(useforplot_vec)) ...
    ', traces were smoothed (0/1): ' num2str(SMOOTHTRACES) '.'])
disp('--------------------------------------------------------')
disp(['bleachrate=' num2str(Bleachrate_effec_mean) ' +- ' num2str(Bleachrate_effec_std) ' 1/sec (mean+-std) (base=2)'])
disp(['half-time T_12=' num2str(T12_effec_mean) ' +- ' num2str(T12_effec_std) ' sec.'])
disp(['offset=' num2str(Offsets_effec_mean) ' +- ' num2str(Offsets_effec_std) ])
disp('--------------------------------------------------------')


%%
% ---------------------------------
% FORMAT IMAGES
% ---------------------------------
set(figure(3),'OuterPosition',[695   499   472   432])
xlim([0 80])
ylim([-0.05 1.2])


set(figure(1),'OuterPosition',[695   499   333   333])
xlim([0 80])
ylim([-0.05 1.2])

set(figure(2),'OuterPosition',[695   499   333   333])
xlim([0 80])
ylim([-0.05 1.2])
%%
%shading image
figure(4)
colormap gray
imagesc(shading)
axis equal
xlim([0 2048]); ylim([0 2048]);
hold on
plot([512 512],[512 1536],'w')
plot([1536 1536],[512 1536],'w')
plot([512 1536],[512 512],'w')
plot([512 1536],[1536 1536],'w')
colorbar
set(figure(4),'OuterPosition',[695   499   356   333])




