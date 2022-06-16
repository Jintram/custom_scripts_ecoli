mydir='D:\ExperimentalDataTodo\ShadingFlatfield\2014-11-11 Flatfield\Flatfield_050ms_bin1\';
d=dir([mydir '*tif*']);
%mydir='D:\ExperimentalDataTodo\Flatfield\Flatfieldprobwrong\';
%d=dir([mydir '*100ms*tif*']);
binning=1;

%get image size
dummyimage=imread([mydir d(1).name]);
imagesize=size(dummyimage);

%stack (3-dim matrix) with 1 flatfield image in each plane
flatfield_complete=zeros(imagesize(1),imagesize(2),length(d));

%load images
for i=1:length(d)
    flatfield_complete(:,:,i)=imread([mydir d(i).name]);
end

myflatfield=uint16(mean(flatfield_complete,3));
myflatfieldFullSize=imresize(myflatfield,binning);
figure(44)
imagesc(myflatfieldFullSize); colorbar;
 %caxis([2*64,2*74])
%caxis([230 260])
colormap(gray)

% name needs to be changed to "flatfield" and then saved together with
% "replace" and "shading"in a matlab-file

%%
 figure(3)
 colormap(gray), 
 imagesc(flatfield)
 colorbar;
 
 flatfieldmean=mean(mean(flatfield))
 newflatfieldmean=mean(mean(myflatfield))