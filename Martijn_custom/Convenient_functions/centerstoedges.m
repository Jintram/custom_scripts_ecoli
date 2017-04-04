function [edges, doubley] = centerstoedges(centers,y)
% equidistant points are assumed

dx=centers(2)-centers(1);
edges = cell2mat(arrayfun(@(x) [centers(x)-dx/2 centers(x)+dx/2], 1:numel(centers),'UniformOutput',0));
doubley = cell2mat(arrayfun(@(x) [y(x) y(x)], 1:numel(y),'UniformOutput',0));


end