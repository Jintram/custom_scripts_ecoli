
% MW 2015/01

% PARAMETER SETTINGS
% ==
Nbins=30;
SMOOTHING=5;

%Look at Dmitry's data

%data(:,1)=randn(100000,1);
%data(:,2)=randn(100000,1);

% Setting up some parameters for the 2d surface stuff
% ===

minx=min(A(:,1));
maxx=max(A(:,1));
miny=min(A(:,2));
maxy=max(A(:,2));

dx=(maxx-minx)/Nbins;
dy=(maxy-miny)/Nbins;

xLinspace = minx:dx:maxx;
yLinspace = miny:dy:maxy;

xCenters = xLinspace(1:end-1)+.5*dx;
yCenters = yLinspace(1:end-1)+.5*dy;

% Plotting original data
figure(1);
plot(A(:,1),A(:,2),'x');

% Making 2d histogram 
% count2d source: http://www.mathworks.com/matlabcentral/fileexchange/29709-function-to-make-a-2d-histogram/content/hist2d.m
% ===
count2d=hist2d(data,xLinspace,yLinspace);
%{
figure(2);
imagesc(count2d); 
%}

histInX = hist(A(:,1),Nbins)
%{
figure(3)
plot(histInX);
%}

[XX,YY] = meshgrid(xCenters,yCenters);
%{
figure(4);
surf(count2d)
%} 

% Normalizing the histogram 
% ===

figure(5);
% Following not used - attempted normalization by histogram of x-direction, see fig 3
myNormalization = meshgrid(histInX,histInX);
normalizedCount2d = count2d./myNormalization;

% Surface plot of 2d-histogram, logarithmic, (+1 to deal w. zero values)
myLogSurface = log(count2d+1)./log(100);
%surf(myLogSurface);
imagesc(xCenters,yCenters,myLogSurface');

% find peaks in y-direction
allpks=zeros(size(myLogSurface)); 
for i=1:size(myLogSurface,1)
    % Pick slice in y-direction
    slice=myLogSurface(i,:);
    % Smooth it
    smoothedslice = smooth(slice,SMOOTHING);
    % Plot it
    %figure(6), hold on, plot(smoothedslice)
    
    % find peaks, put them in allpeaks matrix (holds locations peaks)
    % 1 = peak location
    [currentpks,currentlocs] = findpeaks(smoothedslice);
    currentPeaks = zeros(size(myLogSurface,2),1);
    currentPeaks(currentlocs)=1;
    allpks(i,:) = currentPeaks;
end

% Plot peaks that were found
figure(7)
% surf(allpks)
imagesc(xCenters,yCenters,allpks')

% Watershedding?? 
% ===
% Doesn't work really. Perhaps using y-line localmax peaks as seeds would
% improve it, but that'd take too much time.

negMatrix = -myLogSurface+max(myLogSurface(:));
L = watershed(negMatrix);

figure(10)
subplot(1,2,1)
imshow(negMatrix)
Lrgb = label2rgb(L);
subplot(1,2,2)
imshow(Lrgb)











