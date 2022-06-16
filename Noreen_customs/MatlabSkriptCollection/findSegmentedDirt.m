for i=1:max2(LNsub)
a=regionprops(LNsub==i,'Area');
r=regionprops(LNsub==i,'ConvexArea');
ratio=a.Area/r.ConvexArea;
if ratio<0.85
figure
imagesc(LNsub==i)
title(num2str(i))
end
end