%
imshow2dim=1;
imshowtime=1;
% ****************

mydir='D:\ExperimentalDataTodo\2012-03-08\pos1\';
d_cfp=dir([mydir '*-c-*tif*']);
d_yfp=dir([mydir '*-y-*tif*']);

%get image size
dummyimage=imread([mydir d_cfp(1).name]);
imagesize=size(dummyimage);

%stack (3-dim matrix) with 1 fluor image in each plane
cfp_complete=zeros(imagesize(1),imagesize(2),length(d_cfp));
yfp_complete=zeros(imagesize(1),imagesize(2),length(d_yfp));

%load images
for i=1:length(d_cfp)
    cfp_complete(:,:,i)=imread([mydir d_cfp(i).name]);
end
for i=1:length(d_yfp)
    yfp_complete(:,:,i)=imread([mydir d_yfp(i).name]);
end

cfp_timemean=uint16(mean(cfp_complete,3));
yfp_timemean=uint16(mean(yfp_complete,3));

shading2bin = imresize_old(shading,1/2,'nearest');
      
% **********************************************************************
% 2dim images: spatial
% **********************************************************************
if(imshow2dim==1)

figure(1)
imagesc(cfp_timemean); colorbar;
 %caxis([2*64,2*74])
%caxis([124 135])
colormap(gray)
title('cfp new timemean')


figure(2)
%imagesc(cfp_complete(:,:,1)); colorbar;
imagesc(shading2bin); colorbar;
 %caxis([2*64,2*74])
%caxis([124 135])
colormap(gray)
title('standard shading')


figure(3)
imagesc(yfp_timemean); colorbar;
 %caxis([2*64,2*74])
%caxis([124 135])
colormap(gray)
title('yfp new timemean')


figure(4)
clf
imagesc(double(cfp_timemean)./double(shading2bin)); colorbar;
 %caxis([2*64,2*74])
%caxis([124 135])
colormap(jet)
title('cfp/std_shading','Interpreter','none')

figure(5)
clf
imagesc(double(yfp_timemean)./double(shading2bin)); colorbar;
 %caxis([2*64,2*74])
%caxis([124 135])
colormap(jet)
title('yfp/std_shading','Interpreter','none')

end

% **********************************************************************
% time dependence
% **********************************************************************

if (imshowtime==1)

% average over all pixel values for each frame
cfp_mean_per_frame_dummy=mean(mean(cfp_complete,1),2);
yfp_mean_per_frame_dummy=mean(mean(yfp_complete,1),2);

cfp_mean_per_frame=zeros(size(cfp_mean_per_frame_dummy,3),1);
yfp_mean_per_frame=zeros(size(yfp_mean_per_frame_dummy,3),1);
cfp_mean_per_frame(:)=cfp_mean_per_frame_dummy(1,1,:);
yfp_mean_per_frame(:)=yfp_mean_per_frame_dummy(1,1,:);

figure(6)
plot(cfp_mean_per_frame,'LineWidth',2)
title('average cfp per time')
xlim([0 24])
xlabel('time')
ylim([0 1500])
grid on

figure(7)
plot(yfp_mean_per_frame,'LineWidth',2)
title('average yfp per time')
xlim([0 24])
ylim([0 1500])
xlabel('time')
grid on


end
