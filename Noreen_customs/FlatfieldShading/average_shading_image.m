% ******************************************************************************************
% ******************************************************************************************
% VARIABLES ARE NAMED GFP, BUT SUITABLE ALSO FOR OTHER FLUO CHANNELS

% ------------------------

% *** RAW SHADING IMAGES ******
mydir='D:\ExperimentalDataTodo\ShadingFlatfield\Shading_Micr2\2014-06-23\mCherry_lowconc_20ms_onlytest\';
d_gfp=dir([mydir '*MC*tif*']);
% *** LOAD FINAL FLATFIELD IMAGE (might adjust with  illum time) & REPLACE ******
%load 'D:\SchnitzcellsCurrentVersion\Schnitzcells\fluo_correction_images\Correction_10MHz_GFP_2013_12_18' *flat* replace
% **** WHICH SHADING IMAGES TO SELECT FOR AVERAGING ***
chooseimages=[1:3];

% *** RENAME SHADING AND SAVE MANUALLY WITH FLATFIELD & REPLACE

% ------------------------

%get image size
dummyimage=imread([mydir d_gfp(1).name]);
imagesize=size(dummyimage);

%stack (3-dim matrix) with 1 fluor image in each plane
gfp_complete=zeros(imagesize(1),imagesize(2),length(d_gfp));

for i=1:length(d_gfp)
    gfp_complete(:,:,i)=imread([mydir d_gfp(i).name]);
end


for i=1:size(gfp_complete,3)
    figure(i)
    set(gcf,'WindowStyle','docked')
    clf
    imagesc(gfp_complete(:,:,i))
    colorbar
  % caxis([1000 2200])
end


%%

gfp_choose=gfp_complete(:,:,chooseimages);

gfp_mean=uint16(mean(gfp_choose,3));

%figure(20)
%imagesc(gfp_mean)
%colorbar

gfp_shading_without_flatfield=imresize_old(gfp_mean,2,'nearest');

%gfp_shading_corr_StdFlatfield20111206=gfp_shading_without_flatfield-flatfieldStandard20111206;

%gfp_shading_corr_Flatfield_050ms=gfp_shading_without_flatfield-flatfield050;

%gfp_shading_corr_Flatfield_050ms=gfp_shading_without_flatfield-flatfield;
gfp_shading_corr_Flatfield_new=gfp_shading_without_flatfield-flatfield;

figure(24)
imagesc(gfp_shading_without_flatfield)
title('without flatfield')
 set(gcf,'WindowStyle','docked')
colorbar

%figure(25)
%imagesc(gfp_shading_corr_StdFlatfield20111206)
%title('std flatfield')
%colorbar

figure(26)
imagesc(gfp_shading_corr_Flatfield_new)
title('with flatfield')
 set(gcf,'WindowStyle','docked')
colorbar

%% ******************************************************************************************
% ******************************************************************************************
% ****************OLD ****************
% ******
% COMPARE BETWEEN YFP, CFP and StandardShading Image
load 'D:\ExperimentalDataTodo\2012-04-13shading\collectionShadingDyefiles'  gfp_shading_corr_Flatfield_010ms gfp_shading_corr_Flatfield_25ms cfp_shading_corr_Flatfield_200ms yfp_shading_corr_Flatfield_030ms 
load 'D:\ExperimentalDataTodo\2012-05-08shading\Correction_10MHz_20120508_GFP_NewFlatfield_050ms'  gfp_shading_corr_Flatfield_050ms_2012_05_08
gfp1=double(gfp_shading_corr_Flatfield_25ms); %13.4.
gfp2=double(gfp_shading_corr_Flatfield_010ms); % Maerz
gfp3=double(gfp_shading_corr_Flatfield_050ms_2012_05_08);  %8.5.
cfp1=double(cfp_shading_corr_Flatfield_200ms);
yfp1=double(yfp_shading_corr_Flatfield_030ms);


%normalize
cfp1=cfp1./(mean(mean(cfp1)));
yfp1=yfp1./(mean(mean(yfp1)));
gfp1=gfp1./(mean(mean(gfp1)));
gfp2=gfp2./(mean(mean(gfp2)));
gfp3=gfp3./(mean(mean(gfp3)));



frac_cfp1_yfp1=cfp1./yfp1;
frac_gfp1_gfp2=gfp1./gfp2;
frac_gfp1_gfp3=gfp1./gfp3;
frac_gfp2_gfp3=gfp2./gfp3;
frac_cfp1_gfp1=cfp1./gfp1;
frac_cfp1_gfp2=cfp1./gfp2;
frac_cfp1_gfp3=cfp1./gfp3;
frac_yfp1_gfp1=yfp1./gfp1;
frac_yfp1_gfp2=yfp1./gfp2;
frac_yfp1_gfp3=yfp1./gfp3;

my_caxisrange=[0.9 1.1];

figure(1)
clf
imagesc(frac_gfp1_gfp2)
colorbar
colormap(jet)
%caxis(my_caxisrange)
title('gfp1/gfp2')

figure(2)
clf
imagesc(frac_gfp1_gfp3)
colorbar
colormap(jet)
%caxis(my_caxisrange)
title('gfp1/gfp3')

figure(3)
clf
imagesc(frac_gfp2_gfp3)
colorbar
colormap(jet)
%caxis(my_caxisrange)
title('gfp2/gfp3')


figure(4)
clf
imagesc(gfp1)
colorbar
colormap(jet)
%caxis(my_caxisrange)
title('gfp1')

figure(5)
clf
imagesc(gfp2)
colorbar
colormap(jet)
%caxis(my_caxisrange)
title('gfp2')

figure(6)
clf
imagesc(gfp3)
colorbar
colormap(jet)
%caxis(my_caxisrange)
title('gfp3')



figure(1)
clf
imagesc(frac_yfp1_gfp1)
colorbar
colormap(jet)
caxis(my_caxisrange)
title('yfp1/gfp1')

figure(2)
clf
imagesc(frac_yfp1_gfp2)
colorbar
colormap(jet)
caxis(my_caxisrange)
title('yfp1/gfp2')

figure(3)
clf
imagesc(frac_yfp1_gfp3)
colorbar
colormap(jet)
caxis(my_caxisrange)
title('yfp1/gfp3')


figure(4)
clf
imagesc(frac_cfp1_gfp1)
colorbar
colormap(jet)
caxis(my_caxisrange)
title('cfp1/gfp1')

figure(5)
clf
imagesc(frac_cfp1_gfp2)
colorbar
colormap(jet)
caxis(my_caxisrange)
title('cfp1/gfp2')

figure(6)
clf
imagesc(frac_cfp1_gfp3)
colorbar
colormap(jet)
caxis(my_caxisrange)
title('cfp1/gfp3')




