
%% General
%myim=imread('20140501.tif');         
%myim=imread('20140501_t2.tif');       
myim=imread('20140501_t3.tif');
%myim=imread('dirty_bacteria_t1.tif'); 
%myim=imread('dirty_bacteria_t2.tif'); 

myim=normalize(double(myim)); % always normalize

myiminv=imcomplement(myim);
myimbw=im2bw(myim);

% shape of filter
% ===
%mySE = strel('diamond', 50);
mySE = strel('disk', 11);
%mySE = strel('periodicline', 12,[1 0]);

%% Applying filter - option 1
% TODO THERE IS A BUG SOMEWHERE IN HERE THAT INFLUENCES THE FIXED VARIABLES
% SO RERUN GENERAL SECTION AFTER EXECUTING THIS SECTION
myimin=im2bw(myiminv); % input for this filter

myimero=imerode(myimin, mySE); % erode dirt away
mySE = strel('disk', 11*5);
myimdil=imdilate(myimero,mySE); % dilate bacteria back

figure(1);
subplottight(1,3,1);
imshow(myim,[]);
subplottight(1,3,2);
imshow(myimero,[]);
subplottight(1,3,3);
imshow(myimdil,[]);

%% Tophat
% "substracts opened img from original, enhances contrast"
myimin = myiminv;

mySE = strel('disk', 100);
myimopen = imopen(myimin, mySE); % to illustrate, first step of tophat

% blur on imopen
%h = fspecial('average',100);
%myimopen = imfilter(myimopen, h);

% or to black/white and blur little
myimopen = round(myimopen);
mySE = strel('disk', 25);
myimopen = imerode(myimopen, mySE);
h = fspecial('average',30);
myimopen = imfilter(myimopen, h);


myimtop = myimin-myimopen; % manual tophat (see below for direct matlab)
%myimtop = imtophat(myimin, mySE); % tophat, on original input!

clf;
figure(1);
subplottight(1,3,1)
imshow(myimin,[]);
subplottight(1,3,2)
imshow(myimopen,[]);
subplottight(1,3,3)
imshow(myimtop,[]);

%% skeletonization - leads to a load of artificial branching
% myimin=im2bw(normalize(myimtop),0.2); % tophat as input 
myimin=im2bw(myiminv,0.7);

myimskel=bwmorph(myimin, 'skel', Inf);

clf;
figure(1);
subplottight(1,3,1)
imshow(myim,[]);
subplottight(1,3,2)
imshow(myimin,[]);
subplottight(1,3,3)
imshow(myimskel,[]);


%% Morphological reconstruction (seems not so useful)
myimin = myim; % original is also called "marker "

% Create mask
mask = myimin+.3;

% perform operation
myimrecon = imreconstruct(myimin, mask);

clf;
figure(1);
subplottight(1,3,1)
imshow(myimin,[]);
subplottight(1,3,2)
imshow(mask,[]);
subplottight(1,3,3)
imshow(myimrecon,[]);


%% 3d plot
myimin = myiminv; % original is also called "marker "

% perform operation
myimfill = imfill(myimin);

surf(myimin,'EdgeColor','none','LineStyle','none')%,'FaceLighting','phong');

%% Regional minima
myimin = myiminv;

% Finding max/min
myimmax = imregionalmax(myimin)
myimmin = imregionalmin(myimin)


clf;
figure(1);
subplottight(1,3,1)
imshow(myimin,[]);
subplottight(1,3,2)
imshow(myimmin,[]);
subplottight(1,3,3)
imshow(myimmax,[]);


%% qtdecomp; Some decomposition
% Rather useless for our purposes, but interesting nonetheless.
% (See also example pictures matlab guide.)
close all;
figure()

myims2=myim(1:1024,uint64(((1392-1024)/2)):uint64(((1392-1024)/2+1024))-1);
size(myims2)
myims2=uint8(round(normalize(double(myims2))*256));
S=qtdecomp(uint8(myims2),0.27);
subplottight(1,2,1), imshow(myims2,[])

%%%%%%code from help file
blocks = repmat(uint8(0),size(S));

for dim = [512 256 128 64 32 16 8 4 2 1];    
  numblocks = length(find(S==dim));    
  if (numblocks > 0)        
    values = repmat(uint8(1),[dim dim numblocks]);
    values(2:dim,2:dim,:) = 0;
    blocks = qtsetblk(blocks,S,dim,values);
  end
end

blocks(end,1:end) = 1;
blocks(1:end,end) = 1;

subplottight(1,2,2), imshow(blocks,[])
%%%%%%%%%%%%%%%%%

%% Analyzing texture
% seems useful
figure(1)

myimra = rangefilt(myim);
myimstd = stdfilt(myim);

subplottight(1,3,1), imshow(myim,[])
subplottight(1,3,2), imshow(myimra,[])
subplottight(1,3,3), imshow(myimstd,[])


%% Histogram equalization
% (Note imhist(myim) shows histogram!)
figure(1)

myimeq = histeq(myim,[0:.001:1]);
myimgam = imadjust(myim,[],[],1.5); % gamma=0.5 adjustment

subplottight(1,3,1), imshow(myim,[])
subplottight(1,3,2), imshow(myimeq,[])
subplottight(1,3,3), imshow(myimgam)




%% "Contrast limited histogram equalization"
% seems useful
figure(1)

myimde = decorrstretch(myim);
myimdeli = decorrstretch(myim,'Tol',0.01); % linear contrast stretch also applied

subplottight(1,3,1), imshow(myim,[])
subplottight(1,3,2), imshow(myimde,[])
subplottight(1,3,3), imshow(myimdeli,[])




