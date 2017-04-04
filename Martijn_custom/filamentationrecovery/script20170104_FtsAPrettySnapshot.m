
% Script that makes a snapshot image from the FtsA data set.
% 
% Note: the currently used torquoise color is [48/255, 197/255, 221/255];

%%

load('G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos2crop\segmentation\pos2cropseg137.mat');

% normalize both phsub and yreg
yregnorm = double(yreg); 
yregnorm = (yregnorm-min(yregnorm(:)))./(max(yregnorm(:))-min(yregnorm(:)));
phsubnorm = double(phsub); 
phsubnorm = (phsubnorm-min(phsubnorm(:)))./(max(phsubnorm(:))-min(phsubnorm(:)));

% create contrast-optimalized image for fluor
yregcontrast = double(yreg); 
sortedyreg = sort(yregcontrast(:));
rightPercentileValue = sortedyreg(floor(.999*numel(yregcontrast)));
yregcontrast = (yregcontrast-min(yregcontrast(:)))./(rightPercentileValue-min(yregcontrast(:)));
yregcontrast(yregcontrast>1)=1;
figure(); imshow(yregcontrast,[]);

% create contrast-optimalized image for phsub
phsubcontrast = double(phsub); 
sortedphsub = sort(phsubcontrast(:));
leftPercentileValue = sortedphsub(floor(.01*numel(sortedphsub)));
rightPercentileValue = sortedphsub(floor(.99*numel(sortedphsub)));
phsubcontrast = (phsubcontrast-leftPercentileValue)./(rightPercentileValue-leftPercentileValue);
phsubcontrast(phsubcontrast<0)=0;
phsubcontrast(phsubcontrast>1)=1;
figure(); imshow(phsubcontrast,[]);
figure(); hist(double(phsub(:)),100);
figure(); hist(double(phsubcontrast(:)),100);

%% raw input figures
h1=figure(1); clf; hold on;


subplottight(2,3,1);
imshow(Lc,[]);

subplottight(2,3,2);
imshow(phsub,[]);
subplottight(2,3,3);
imshow(yreg,[]);

subplottight(2,3,5);
imshow(phsubcontrast,[]);
subplottight(2,3,6);
imshow(yregcontrast,[]);


%% Get perimeter (outlines) of cells

%Code stolen from
% edit PN_imshowlabel

L=Lc;

cellBodies=zeros(size(L));

for cellIdx = 1:max(max(L))

    workImg=zeros(size(L));
    workImg(L==cellIdx)=1;

    cellBodies = cellBodies+imerode(workImg,strel('disk',1,4));
end

perimImg(L>0)=1;
perimImg(cellBodies>0) = 0;

% Put perimeters in phase image
nonZeroIdx = perimImg>0; %perimImg(:,:,1)>0; | perimImg(:,:,2)>0 | perimImg(:,:,3)>0; % some RGB values have a 0 in them, e.g. red = [1,0,0]
%nonZeroIdx = [nonZeroIdx,nonZeroIdx,nonZeroIdx];

%% yreg w. outlines

h2=figure(2); clf; hold on;

% initialize image
outim = double(yreg); 
% normalize image
outim = (outim-min(outim(:)))./(max(outim(:))-min(outim(:)));
% introduce outlines (perimiters)
outim(perimImg>0)=.5;
% fade out non-cell background
outim(Lc==0)=outim(Lc==0)*.5;

% 
imshow(outim,[]);

%% combination of phase and fluor

h4=figure(4); clf; hold on;

% Make 3d phaseimg
FTSABLENDRATIO = .5;
outim2=zeros([size(phsubcontrast),3]); % phase as base          
outim2(:,:,1) = phsubcontrast.*(1-FTSABLENDRATIO)+yregcontrast.*255/255.*FTSABLENDRATIO;
outim2(:,:,2) = phsubcontrast.*(1-FTSABLENDRATIO)+yregcontrast.*255/255.*FTSABLENDRATIO;
outim2(:,:,3) = phsubcontrast.*(1-FTSABLENDRATIO)+yregcontrast.*0/255.*FTSABLENDRATIO; 

% introduce outlines (perimiters)
%{
zeroesImgSize=zeros(numel(yregnorm),1);
ThePerimIndx = perimImg(:)>0;
outim2(find([ThePerimIndx; zeroesImgSize; zeroesImgSize]'))=48/255; % 48/255;%240/255;
outim2(find([zeroesImgSize; ThePerimIndx; zeroesImgSize]'))=197/255; %197/255;%203/255;
outim2(find([zeroesImgSize; zeroesImgSize; ThePerimIndx]'))=221/255;  %221/255;%34/255;
%}

%outim2 = phsubnorm.*.5+yregnorm.*.5;

%
imshow(outim2,[]);

%% yregcontrast with coloured cell outlines

h5=figure(5); clf; hold on;

% Make 3d phaseimg
outim=zeros([size(yregcontrast),3]); % phase as base          
outim(:,:,1) = yregcontrast;
outim(:,:,2) = yregcontrast;
outim(:,:,3) = yregcontrast;

% introduce outlines (perimiters)
zeroesImgSize=zeros(numel(yregcontrast),1);
ThePerimIndx = perimImg(:)>0;
outim(find([ThePerimIndx; zeroesImgSize; zeroesImgSize]'))=48/255; % 48/255;%240/255;
outim(find([zeroesImgSize; ThePerimIndx; zeroesImgSize]'))=197/255; %197/255;%203/255;
outim(find([zeroesImgSize; zeroesImgSize; ThePerimIndx]'))=221/255;  %221/255;%34/255;

% fade out non-cell background
TheFadeoutIndx = Lc(:)==0;
outim(find([TheFadeoutIndx; zeroesImgSize; zeroesImgSize]'))=...
    .5.*outim(find([TheFadeoutIndx; zeroesImgSize; zeroesImgSize]'));
outim(find([zeroesImgSize; TheFadeoutIndx; zeroesImgSize]'))=...
    .5.*outim(find([zeroesImgSize; TheFadeoutIndx; zeroesImgSize]'));
outim(find([zeroesImgSize; zeroesImgSize; TheFadeoutIndx]'))=...
    .5.*outim(find([zeroesImgSize; zeroesImgSize; TheFadeoutIndx]'));


%outim(Lc==0)=outim(Lc==0)*.5;

% 
imshow(outim,[]);
hFig3a = h5;


