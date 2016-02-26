function h = subplottightmargin(n,m,i,margin)
% function h = subplottightmargin(n,m,i,margin)
% MW: function that produces tight subplot    
% Source: 
% http://www.briandalessandro.com/blog/how-to-make-a-borderless-subplot-of-images-in-matlab/

[c,r] = ind2sub([m n], i);

marginc = 1/m*margin;
marginr = 1/n*margin;

ax = subplot('Position', [ (c-1)/m, 1-(r)/n, ...
                           1/m-marginc, 1/n-marginr]);

if(nargout > 0)
  h = ax;
end
   

end