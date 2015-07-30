
function mymap = loadcolormapdjk(L)
% Based on PN_imshowlabel.m; adapted by MW 29-8-2014
%
% Create DJK's color map
%
% L should contain matrix with segm. bacteria
%
% function mymap = loadcolormapdjk(L)

L2 = mod(L,255)+2;
L2(L==0) = 1;
% M is the maximum color table entry, at most 256 colors
M = min(max2(L)+2,256);

% create a color map
mymap = DJK_hsv(M); % DJK 071207
% explicitly set the colormap's first entry to black for background
mymap(1,:)=[0 0 0];
%if p_internal.randomize
  % get sequence of random integers in range [1,maxcolors-1]
  [s,I] = sort(rand(M-1,1));  
  % randomly reorder mymap color entries [2,maxcolors]
  mymap(2:end,:) = mymap(I+1,:);
%end
mymap = [mymap ; 1 1 1]; %add white

end