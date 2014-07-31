% create a set of fluo image with half the size on second diagonal (Raute) and a grey background to fill the full
% size -> check if rescaling is correctly inserted into every part of schnitzcells
% rescentering skipped because figs should be in centre anyways
%METADATA HAS TO BE ADDED MANUALLY!

%*****************ADJUST************************
rescalefactor=[1 0.5];
mydir='D:\ExperimentalDataTodo\2012-01-27\pos5crop\images\';
d=dir([mydir, '*-c-*']);
mysavedir='D:\ExperimentalDataTodo\2012-01-27\pos5crop\images\extraCFPscale\'; % make sure it exists!
%for debugging: show figures
showfig=1;
    
%***********************************************
for runfluor=1:length(d)
    origimage=imread([mydir, d(runfluor).name]);
    originfo= imfinfo([mydir, d(runfluor).name]);
    
    disp(['loaded ', d(runfluor).name])




if length(rescalefactor)==1
    disp(['Only one rescale factor given instead of two. Will perform homogeneous stretching along both diagonals']);
    rescalefactor=[rescalefactor, rescalefactor];
end
% stretch factors
s1=rescalefactor(1);
s2=rescalefactor(2);


% input image
%origimage=checkerboard(10,2); debugging
%origimage=imread('D:\DummyExp\2012-fl-uo\TestData\pos5crop-c-035.tif');
if showfig==1
    figure(1);
    clf
    colormap('gray');
    %imshow(origimage,'Border','loose')
    imagesc(origimage)
    hold on
    axis on
end

%max coordinates
dims=size(origimage); 
d1=dims(1); d2=dims(2);
%diagonal vectors from start to end
% v1=top left -> bottom right. v2=bottom left -> top right
v1=[d2-1,d1-1]; %(1=initial coordinate in upper left corner (not =0!))
v2=[d2-1,1-d1];
% better readability (like in script)
x1=v1(1);y1=v1(2);x2=v2(1);y2=v2(2); %y2<0. x2>0

% transformation matrix + elements;
c=x1*y2-x2*y1; % pre constant
t11=s1*x1*y2-s2*x2*y1;
t12=y1*y2*(s1-s2);
t21=x1*x2*(s2-s1);
t22=s2*x1*y2-s1*x2*y1;
T= 1/c*[t11  t12;
    t21  t22;
    0     0]; % no shift

t_aff = maketform('affine',T);
image_affine = imtransform(origimage,t_aff,'FillValues',0.3,'XYScale',1);

if showfig==1
    figure(2);
    hold off
    clf
    colormap('gray');
    %imshow(image_affine,'Border','loose')
    imagesc(image_affine)
    hold on
    axis on
    %axis image
end


% save rescaled fluor image
% *****************************************************************************
% AT THE MOMENT ARBITRARY FILE NAME
% *****************************************************************************
% at the moment centerfluor is unint16, check later!
%image_affine = imresize_old(image_affine,0.5,'nearest'); % not necessary,
% % since not binning corected
savename=[mysavedir d(runfluor).name];
imwrite(image_affine,savename ,'TIFF','Compression','none','Description',originfo.ImageDescription)
  disp(['saved ', d(runfluor).name])

end
