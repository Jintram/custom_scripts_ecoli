% checks fluctuations in lamp intensity over time by directly loading all fluo
% images. then flatfield correction. then shading correction
% different averages of the brightest pixels are calculated (e.g. 1st till
% 1000th brightest etc). plot over time.
% does not use any segmentation put only the pure images

clear mydir dir_fluoimages numim dummyimage imagesize fluoimages_complete flatfield


% ***********
% ADJUST!!!
% if # of observed beads changes -> restrict used framenumbers. otherwise,
% suddenly e.g. 4 oinstead of 3 bright bead-centres are in the field of
% view -> increase in average intensity!!)
%frameRange=[0, 778]; %min and max allowed frame number (pos6. 2013-04-16)
%frameRange=[0 593]; %min and max allowed frame number (pos7. 2013-04-16)
%mydir='D:\ExperimentalDataTodo\2013-04-16\pos6\images\';  %careful with folders!
%mysaveDir='D:\ExperimentalDataTodo\2013-04-16\pos6crop\analysis\schnitzcells\';
%mydir='D:\ExperimentalDataTodo\2013-04-16\pos7\images\';  %careful with folders!
%mysaveDir='D:\ExperimentalDataTodo\2013-04-16\pos7crop\analysis\schnitzcells\';
mydir='D:\ExperimentalDataTodo\2013-04-26\focus1\2013-04-26\pos8\';
mysaveDir='D:\ExperimentalDataTodo\2013-04-26\focus1\2013-04-26\pos8\';
frameRange=[0 8];
% ***********


% ---------------------------------------------------------
% LOAD IMAGES
% ---------------------------------------------------------

dir_fluoimages=dir([mydir '*-r-*tif*']);  %adjust fluo color name!
numim=length(dir_fluoimages);
% assumes binning of=1

% LOAD A FLATFIELD IMAGE
load 'D:\SchnitzcellsCurrentVersion\Schnitzcells\fluo_correction_images\Correction_10MHz_20130418_GFP_flat100msDigi20Mhz' flatfield shading
flatfield=double(flatfield);
shading=double(shading);

%get image size
dummyimage=imread([mydir dir_fluoimages(1).name]);
imagesize=size(dummyimage);

%stack (3-dim matrix) with 1 fluor image in each plane
fluoimages_complete=zeros(imagesize(1),imagesize(2),numim);

for i=1:numim
    fluoimages_complete(:,:,i)=imread([mydir dir_fluoimages(i).name]);
    %subtract flatfield
    fluoimages_complete(:,:,i)=fluoimages_complete(:,:,i)-flatfield;
    fluoimages_complete(:,:,i)=fluoimages_complete(:,:,i)./shading*mean2(shading);
end


%get frame numbers of images (to label x-axis and to compare with
%schnitzcells-traces)
framenumbers=[];
for i=1:length(dir_fluoimages)
    curname=dir_fluoimages(i).name;
    %ADJUST NAME IF COLOR CHANGED!!!!!
    start=strfind(curname,'-r-');
    strframe=curname(start+3:start+5);
    framenumbers=[framenumbers; str2num(strframe)];
end


%select allowed framerange
% if optimization necessary: do this before image loading.
frameidx=find(framenumbers>=frameRange(1) & framenumbers<=frameRange(2));
fluoimages_complete=fluoimages_complete(:,:,frameidx);
%recalculate #images and relevant framenumbers
numim=length(frameidx);
framenumbers=framenumbers(frameidx);

%for i=1:size(fluorimages_complete,3)
%    figure(i)
%    clf
%    imagesc(fluorimages_complete(:,:,i))
%    colorbar
%  % caxis([1000 2200])
%end



% ---------------------------------------------------------
% CALCULATE AVERAGE PIXEL INTENSITY TRACES 
% ---------------------------------------------------------

% calculate several pixel averages
mean1to10000=zeros(numim,1); %stores for each fluo image the mean of the 1st till 10000th brighest pixel
mean1000to10000=zeros(numim,1); % in case of shot noise in brightest pixels
mean1to2000=zeros(numim,1);
mean1000to2000=zeros(numim,1);
meanMin100=zeros(numim,1); % mean of all pixels with minimum value 100 (after flatfiedl corr!)
sumMin100=zeros(numim,1); % sum of all pixels with minimum value 100 (after flatfiedl corr!)

clear pixelssorted
% loop over all frames
for fr=1:numim
    currentimage=fluoimages_complete(:,:,fr);
    pixelssorted=sort(currentimage(:));
    pixelssorted=pixelssorted(end:-1:1);
    mean1to10000(fr)=mean(pixelssorted(1:10000));
    mean1000to10000(fr)=mean(pixelssorted(1000:10000));
    mean1to2000(fr)=mean(pixelssorted(1:2000));
    mean1000to2000(fr)=mean(pixelssorted(1000:2000));
    brightidx=find(pixelssorted>=100);
    brightpixels=pixelssorted(brightidx);
    meanMin100(fr)=mean(brightpixels);
    sumMin100(fr)=sum(brightpixels);
end
clear currentimage pixelssorted


% ---------------------------------------------------------
% PLOTTING
% ---------------------------------------------------------

% automatically dock figures
set(0,'DefaultFigureWindowStyle','docked') 

%histogram of pixel intensity for x'th image
testnr=2;
testimage=fluoimages_complete(:,:,testnr);
pixelssorted=sort(testimage(:));
pixelssorted=pixelssorted(end:-1:1);
figure(1)
clf
subplot(2,1,1)
plot(pixelssorted)
title(['pixel distribution of image ', num2str(testnr)])
xlabel('x''th pixel')
ylabel('its value')
grid on
subplot(2,1,2)
plot(pixelssorted)
hold on
plot(ones(500,1)*1000,1:1:500,'--k')
plot(ones(500,1)*2000,1:1:500,'--k')
xlim([0 15000])
title('subplot')
xlabel('x''th pixel')
ylabel('its value')
grid on

clear testimage pixelssorted



% plot averages over time
%figure(2)
%clf
%plot(framenumbers,mean1to10000,'.-')
%title('mean1to10000')

%figure(3)
%clf
%plot(framenumbers,mean1000to10000,'.-')
%title('mean1000to10000')

%figure(4)
%clf
%plot(framenumbers,mean1to2000,'.-')
%title('mean1to2000')

%figure(5)
%clf
%plot(framenumbers,mean1000to2000,'.-')
%title('mean1000to2000')

% combined plot of different pixel averaging methods
figure(5)
clf
plot(framenumbers,mean1to10000/mean(mean1to10000),'.-b')
hold on
plot(framenumbers,mean1000to10000/mean(mean1000to10000),'.-c')
plot(framenumbers,mean1to2000/mean(mean1to2000),'.-r')
plot(framenumbers,mean1000to2000/mean(mean1000to2000),'.-m')
grid on
legend('mean1to10000','mean1000to10000','mean1to2000','mean1000to2000')
title('combined & normalized')
ylabel('pixel intensity')
xlabel('frame number')

%save image
figureFileName=['TimeTraces_XYFrames_vs_PixelIntensity_fixPxNr'];
saveSameSize(gcf,'file',[mysaveDir figureFileName '.png'], 'format', 'png');

% combined plot of different pixel averaging methods
figure(6)
clf
plot(framenumbers,meanMin100/mean(meanMin100),'.-k')
hold on
plot(framenumbers,sumMin100/mean(sumMin100),'.-','Color', [1 0.5 0])
grid on
legend('meanMin100','sumMin100')
title('all brightest pixels')
ylabel('pixel intensity')
xlabel('frame number')

%save image
figureFileName=['TimeTraces_XYFrames_vs_PixelIntensity_allBrightPx'];
saveSameSize(gcf,'file',[mysaveDir figureFileName '.png'], 'format', 'png');

%set(0,'DefaultFigureWindowStyle','normal') 

%%

% ---------------------------------------------------------
% CALCULATE TRACES FOR EACH BRIGHT SPOT INDIVIDUALLY
% 1st STEP DETERMINE REGIONS OF INTEREST
% (no segmentation existent!)
% ---------------------------------------------------------

% load first image and a late image and define area of interest (coordinates) for each bead
% (i.e. only one bead is contained in that area)
% find latest image (pad might move!) for which still the same rectangular
% regions of interest can be defined

framenumberofinterest=5;

figbla=figure(10);
clf
subplot(2,1,1)
imagesc(max2(fluoimages_complete(:,:,1))-fluoimages_complete(:,:,1))
title('first image')
grid on
xlabel('x (1st & 2nd in regcoord)')
ylabel('y (3rd & 4th in regcoord)')
colormap ('gray')

subplot(2,1,2)
idxframe=find(framenumbers==framenumberofinterest);
if ~isempty(idxframe)
    imagesc(max2(fluoimages_complete(:,:,idxframe))-fluoimages_complete(:,:,idxframe))
    title(['frame ' num2str(framenumberofinterest) ' image (idx:' num2str(idxframe) ')'])
    grid on
end

% axis structure
% 0 -------------> (x)
% |
% |
% |
% (y)

% ******** ADJUST ******************
% frame rane for which cst region of interest rectangles for each bead can
% be defined (pad must not move too much!)
frameRangeSingle=[0 593];
% REGIONS OF INTEREST
regcoord=[]; % lower(x) - upper(x) - lower(y) - upper(y) (1st bead)
                       % lower(x) - upper(x) - lower(y) - upper(y) (2nd bead)
%regcoord= [150 400 150 400;
%          200 500 500 800;
%          700 900 600 900];
regcoord= [900 1000 300 450; 600 800 400 600         ];
          


% ******** ENDADJUST ***************
numberofregions=size(regcoord,1);


%select allowed framerange
frameSingleidx=find(framenumbers>=frameRangeSingle(1) & framenumbers<=frameRangeSingle(2));

figure(figbla)
subplot(2,1,1)
hold on
for i=1:numberofregions
    rectangle('Position',[regcoord(i,1),regcoord(i,3),regcoord(i,2)- ...
        regcoord(i,1), regcoord(i,4)-regcoord(i,3)])
end
subplot(2,1,2)
hold on
for i=1:numberofregions
    rectangle('Position',[regcoord(i,1),regcoord(i,3),regcoord(i,2)- ...
        regcoord(i,1), regcoord(i,4)-regcoord(i,3)])
end


%%
% ---------------------------------------------------------
% 2nd STEP: DETERMINE TIME TRACES OF INDIVIDUAL SPOTS
% ---------------------------------------------------------

allROIDatamean1to2000=zeros(0,3); % mean of 1st to 2000th brighest pixel in ROI
allROIDatameanMin100=zeros(0,3); % mean of all pixels in ROI with >100 intensity
%allROIDataSumPixels=zeros(0,3);  % sums all pixels in ROI
    % THE SUM OF ALL PIXELS TURNS OUT TO BE A BAD MEASURE: super sensitive
    % to size of ROI (=# background pixels)
    
% contains complete x&y data of all ROIs (region of interests)
% ROI1  ---  x (framenr1) --- mean1to2000_1
% ROI1  ---  x (framenr2) --- mean1to2000_2
% ...
% ROI2  ---  x (framenr2) --- mean1to2000_1
% ROI2  ---  x (framenr2) --- mean1to2000_2
% etc etc

% loop over all ROIs
for roirun=1:numberofregions
    subimages=fluoimages_complete(regcoord(roirun,3):regcoord(roirun,4),regcoord(roirun,1):regcoord(roirun,2),frameSingleidx);
    for fr=1:length(frameSingleidx)
        curframe=framenumbers(frameSingleidx(fr));
        clear mean1to2000single meanMin100single SumPixelssingle
        currentimage=subimages(:,:,fr);
        pixelssorted=sort(currentimage(:));
        pixelssorted=pixelssorted(end:-1:1);
        % get time trace data point
        mean1to2000single=mean(pixelssorted(1:2000));
        brightidx=find(pixelssorted>=100);
        brightpixels=pixelssorted(brightidx);
        meanMin100single=mean(brightpixels);
        %SumPixelssingle=sum(sum(currentimage));
        %add to matrix
        allROIDatamean1to2000=[ allROIDatamean1to2000; roirun, curframe, mean1to2000single];
        allROIDatameanMin100=[ allROIDatameanMin100; roirun, curframe, meanMin100single];
        %allROIDataSumPixels= [ allROIDataSumPixels; roirun, curframe, SumPixelssingle];
    end
end
clear currentimage pixelssorted curframe fr mean1to2000 meanMin100 SumPixels

%%
% ---------------------------------------------------------
% PLOT INDIVIDUAL TRACES OVER TIME AND THEIR AVERAGE
% ---------------------------------------------------------

fig2=figure(21);
set(gcf,'WindowStyle','docked')
clf
hold on
xlabel('frame number')
title('individual ROI traces. normalized to 1')
grid on

%define colors
myColormap=[1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0]; 
number_of_colors=size(myColormap,1);

% plot timetraces for all ROIs
% intensity measure (currently): meanMin100
ylabel('meanMin100','Interpreter','None')

for roirun=1:numberofregions
    idxroirun=find(allROIDatameanMin100(:,1)==roirun);
    currentframes=allROIDatameanMin100(idxroirun,2);
    currentyvalues=allROIDatameanMin100(idxroirun,3);
    %determine color
    coloridx=mod((roirun-1),(number_of_colors-1))+1;
    %plot
    plot(currentframes,currentyvalues/mean(currentyvalues),'.-','Color',myColormap(coloridx,:))
    legend('1','2')
end

%save image
figureFileName=['TimeTraces_XYFrames_vs_PixelIntensity_IndivROITraces'];
saveSameSize(gcf,'file',[mysaveDir figureFileName '.png'], 'format', 'png');




%%
% ---------------------------
% QUICKCHECK
% ---------------------------

a1=allROIDatamean1to2000(1:74,3);
a2=allROIDatamean1to2000(75:148,3);
a3=allROIDatamean1to2000(149:end,3);
b1=allROIDatameanMin100(1:74,3);
b2=allROIDatameanMin100(75:148,3);
b3=allROIDatameanMin100(149:end,3);
c1=allROIDataSumPixels(1:74,3);
c2=allROIDataSumPixels(75:148,3);
c3=allROIDataSumPixels(149:end,3);

corr(a1,c1)
corr(a2,c2)
corr(a3,c3)

figure
plot(a1/mean(a1),'.-b')
hold on
plot(b1/mean(b1),'.-r')
plot(c1/mean(c1),'.-g')
title('first ROI')

figure
plot(a2/mean(a2),'.-b')
hold on
plot(b2/mean(b2),'.-r')
plot(c2/mean(c2),'.-g')
title('second ROI')

figure
plot(a3/mean(a3),'.-b')
hold on
plot(b3/mean(b3),'.-r')
plot(c3/mean(c3),'.-g')
title('third ROI')

figure
plot(a1/mean(a1),'.-b')
hold on
plot(a2/mean(a2),'.-r')
plot(a3/mean(a3),'.-g')
title('mean1to2000')

figure
plot(b1/mean(b1),'.-b')
hold on
plot(b2/mean(b2),'.-r')
plot(b3/mean(b3),'.-g')
title('meanMin100')


figure
plot(c1/mean(c1),'.-b')
hold on
plot(c2/mean(c2),'.-r')
plot(c3/mean(c3),'.-g')
title('sumpixels')
