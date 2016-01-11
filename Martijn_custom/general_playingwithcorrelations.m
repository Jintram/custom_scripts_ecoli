% Simple script to get some feel for correlation funtions
% ===
%
% This script calculate corr or cov values for a cloud.
%
% Note that Noreen has sent me two websites with two nice explanations when
% the slope indeed is the same as the cov, and when not.
%
% - http://stats.stackexchange.com/questions/22718/what-is-the-difference-between-linear-regression-on-y-with-x-and-x-with-y/22721#22721
% - http://stats.stackexchange.com/questions/32464/how-does-the-correlation-coefficient-differ-from-regression-slope
%
% (See also word file where I copied the contents of those sites.)

fi = 1; % figure index

%% Simple linear relationsship
% ===
x = [0:.01:1];
y = 5*x+2;
figure(fi); clf; plot(x,y);
Rlin = corr(x',y')
%corrcoef(x',y')
fi = fi+1;

%% The outcome of this is one, which can be derived easily,
% see notes 19-04-2015 (p. 2). A manual calculation of the R is as follows:
% ===
% Note:
% corr(x(t),y(t)) = 
varXY = 0; varX = 0; varY = 0;
xNorm = x-mean(x);
yNorm = y-mean(y);
for t=1:numel(x)
    varXY = varXY + xNorm(t)*yNorm(t);    
    varX = varX+xNorm(t)^2;
    varY = varY+yNorm(t)^2;
end
RXY = varXY / sqrt(varX*varY)

% Note that you can calculate the correlation coefficients easily:
betaY = varXY / varX
betaX = varXY / varY
% Note furthermore that these relations are derived by wanting to find the
% minimum of the sum of squared errors; this can be done by looking for the
% derivative being zero. Above relations then follow. See also:
% http://en.wikipedia.org/wiki/Ordinary_least_squares
% http://en.wikipedia.org/wiki/Proofs_involving_ordinary_least_squares#Least_squares_estimator_for_.CE.B2
% (And this was also mentioned in my programming course - see lecture notes 
% thereof.)

%% What happens if you have a non-linear relation?
% E.g. a parabola
% ===
x = [-1:.01:1];
y = -1*x.^2;
figure(fi); clf; plot(x,y);
Rpara = corr(x',y')
fi = fi+1;

% % Rpara =
% %  3.1614e-18
% So apparently the correlation disappears completely for a parabola.
% Which makes sense since there is no straight line describing this 
% relationship.

%% Another example
% ===
x = [-1:.01:1];
y = 4*x.^3+x.^4;
figure(fi); clf; plot(x,y);
Rarbi = corr(x',y')
fi = fi+1;












