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
EXPORTFOLDER='\\storage01\data\AMOLF\users\wehrens\Latex3\Thesis\Chapter2_Methods\Figures\MatlabExport\';

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
N = numel(x);
for tf=1:N
    % sum
    varXY = (varXY + xNorm(tf)*yNorm(tf));
    varX  = (varX+xNorm(tf)^2);
    varY  = (varY+yNorm(tf)^2);
end
% normalize
varXY = varXY/N; varX = varX/N; varY = varY/N;
% correlation coefficient
RXY = varXY / sqrt(varX*varY)

% Note that you can calculate the least square fit coefficients easily:
betaY = varXY / varX
betaX = varXY / varY

bSanityY = Rlin * sqrt(var(x)*var(y)) / var(x)
%
% For y=betay*x+a1 and x=betax*y+a2
% 
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


%% An example with a cloud and some noise, linear relation

x = [-1:.01:1];
%y = 4*x.^3+x.^4;
y = 5*x+2;
mysigmax=1;
mysigmay=1;
xrand=[];
yrand=[];
for n=1:250
    idx=int16(rand()*(numel(x)-1)+1);
    xrand(end+1)=x(idx)+normrnd(0,mysigmax);
    yrand(end+1)=y(idx)+normrnd(0,mysigmay);
end

hEx1=figure(2); clf; hold on;
scatter(xrand,yrand,'LineWidth',2);
plot(x,y,'-k','LineWidth',2);

%legend('Added noise','Given')

xlabel('f');
ylabel('g');

[Rnoisy1,pnoisy1] = corr(xrand',yrand')
text(-4,9+2,['R = ' sprintf('%0.3f', Rnoisy1)]);
text(-4,7+2,['p = ' num2str(pnoisy1)]);
%text(-4,7+2,['p = ' sprintf('%0.7f', pnoisy1)]);



% Note that calculating b doesn't work if the input is noisy (try e.g.
% sigma=.01, then it works well).
betaY = Rnoisy1 * sqrt(var(xrand)*var(yrand)) / var(xrand)
bSanityY = corr(xrand',yrand') * sqrt(var(xrand)*var(yrand)) / var(xrand)
% lsqlin(yrand,2,xrand,100) % lsqlin(C,d,A,b) % Not sure what input I
% should give here..

% Now also add the mutual information
[pdf1, pdf2, pdf1pdf2, pdf12, dx, dy] = getpfds12MW([xrand;yrand]',256,1);
mutualInformation = entropyMW(pdf1)+entropyMW(pdf2)-entropyMW(pdf12)
figure(2);
text(.5,-6,['I = ' sprintf('%0.3f', mutualInformation)]);


%MW_makeplotlookbetter(14,[],[6,4]);
MW_makeplotlookbetter(14,[],[(12.8-.4*3)/3,4]);

ylim([-10,10]+2);
xlim([-5,5]);


%% An example with a cloud and some noise, linear relation

x = [-1:.01:1];
y = -7*(x).^2+4;%+x.^4;
%y = 5*x+2;
mysigmax=.25;
mysigmay=.25;
xrand=[];
yrand=[];
idxregister=[];
for n=1:250
    idx=int16(rand()*(numel(x)-1)+1);
    xrand(end+1)=x(idx)+normrnd(0,mysigmax);
    yrand(end+1)=y(idx)+normrnd(0,mysigmay);
    idxregister(end+1)=idx;
end
% consistency check:
% figure; hist(idxregister,100)

hEx2=figure(3); clf; hold on;

scatter(xrand,yrand,'LineWidth',2);

[bandwidth,density,X,Y] = kde2d([xrand', yrand']); 
[C, l1] = contour(X,Y,density,5,'-','LineWidth',2,'Color',[.5 .5 .5]);

plot(x,y,'-k','LineWidth',2);

%legend('Added noise','Given')

xlabel('f');
ylabel('g');

[Rnoisy1,pnoisy1] = corr(xrand',yrand')
text(-2.5,9+2,['R = ' sprintf('%0.3f', Rnoisy1)]);
text(-2.5,7+2,['p = ' num2str(pnoisy1)]);
%text(-4,7+2,['p = ' sprintf('%0.7f', pnoisy1)]);

% Now also add the mutual information
[pdf1, pdf2, pdf1pdf2, pdf12, dx, dy] = getpfds12MW([xrand;yrand]',256,1);
mutualInformation = entropyMW(pdf1)+entropyMW(pdf2)-entropyMW(pdf12)
figure(3);
text(.5,-6,['I = ' sprintf('%0.3f', mutualInformation)]);

%MW_makeplotlookbetter(14,[],[6,4]);
MW_makeplotlookbetter(14,[],[(12.8-.4*3)/3,4]);

ylim([-10,10]+2);
xlim([-3,3]);

%% Same as previous, but sigma=1

x = [-1:.01:1];
y = -7*(x).^2+4;%+x.^4;
%y = 5*x+2;
mysigmax=1;
mysigmay=1;
xrand=[];
yrand=[];
idxregister=[];
for n=1:250
    idx=int16(rand()*(numel(x)-1)+1);
    xrand(end+1)=x(idx)+normrnd(0,mysigmax);
    yrand(end+1)=y(idx)+normrnd(0,mysigmay);
    idxregister(end+1)=idx;
end
% consistency check:
% figure; hist(idxregister,100)

hEx3=figure(4); clf; hold on;

scatter(xrand,yrand,'LineWidth',2);

[bandwidth,density,X,Y] = kde2d([xrand', yrand']); 
[C, l1] = contour(X,Y,density,5,'-','LineWidth',2,'Color',[.5 .5 .5]);

plot(x,y,'-k','LineWidth',2);

%legend('Added noise','Given')

xlabel('f');
ylabel('g');

[Rnoisy1,pnoisy1] = corr(xrand',yrand')
text(-2.5,9+2,['R = ' sprintf('%0.3f', Rnoisy1)]);
text(-2.5,7+2,['p = ' num2str(pnoisy1)]);
%text(-4,7+2,['p = ' sprintf('%0.7f', pnoisy1)]);

% Now also add the mutual information
[pdf1, pdf2, pdf1pdf2, pdf12, dx, dy] = getpfds12MW([xrand;yrand]',256,1);
mutualInformation = entropyMW(pdf1)+entropyMW(pdf2)-entropyMW(pdf12)
figure(4);
text(.5,-6,['I = ' sprintf('%0.3f', mutualInformation)]);

%MW_makeplotlookbetter(14,[],[6,4]);
MW_makeplotlookbetter(14,[],[(12.8-.4*3)/3,4]);

ylim([-10,10]+2);
xlim([-3,3]);

% Another consistensy check
% figure; histogram(normrnd(0,mysigmax,250,1),30)



%% 

if exist('SAVEPLEASE','var')
    figure(hEx1); MW_makeplotlookbetter(14,[],[3.8,4]);
    saveas(hEx1,[EXPORTFOLDER 'SVG_technical_scatter1.svg']);
    saveas(hEx1,[EXPORTFOLDER 'TIF_technical_scatter1.tif']);
    saveas(hEx1,[EXPORTFOLDER 'FIG_technical_scatter1.fig']);

    figure(hEx2); MW_makeplotlookbetter(14,[],[3.8,4]);
    saveas(hEx2,[EXPORTFOLDER 'SVG_technical_scatter2.svg']);
    saveas(hEx2,[EXPORTFOLDER 'TIF_technical_scatter2.tif']);
    saveas(hEx2,[EXPORTFOLDER 'FIG_technical_scatter2.fig']);
    
    figure(hEx3); MW_makeplotlookbetter(14,[],[3.8,4]);
    saveas(hEx3,[EXPORTFOLDER 'SVG_technical_scatter3.svg']);
    saveas(hEx3,[EXPORTFOLDER 'TIF_technical_scatter3.tif']);
    saveas(hEx3,[EXPORTFOLDER 'FIG_technical_scatter3.fig']);
end

%% a delayed relationship between two diffusive processes, illustrated 
% by a random walk process and a process that is part random walk and part
% related to the first walk.

mylinecolors=linspecer(4);
DELTA=90; EXTRA=DELTA+50;

n=1:1000; 
nprime=100:1000;

% simple random walk processes
f       =cumsum(arrayfun(@(i) rand()-.5,n)); % f is defined for all n
gpart1  =cumsum(arrayfun(@(i) rand()-.5,nprime)); % g is only defined for nprime

% now create g, as sum of gpart1 and gpart2, still only defined at nprime,
% nprime(1) matchines g(1), and nprime(end) matching g(end).
gpart2=f(nprime-DELTA).^2
g=.8*gpart1+.2*gpart2;

% plot signals
% hdelaysignal=figure(5); clf; hold on;
hdelay=figure(5); clf;
subplot(1,3,1); hold on;
plot(nprime-100,f(nprime),'Color',mylinecolors(1,:),'LineWidth',2)
plot(nprime-100,g,'Color',mylinecolors(2,:),'LineWidth',2)

% cosmetic
xlabel('t');
ylabel('Function value');
MW_makeplotlookbetter(14,[],[12.8,4]);

% plot scatter 1
%hdelayscatter=figure(6); clf; hold on;
subplot(1,3,2); hold on;
plot(f(nprime),g,'.','LineWidth',2,'Color',mylinecolors(3,:),'MarkerSize',3^2);
%plot(f(n+DELTA),g(n),'.','LineWidth',2,'Color',mylinecolors(2,:),'MarkerSize',3^2);

% cosmetic
%legend('f(n) vs g(n)','f(n) vs g(n+\tau)')
xlabel('f(t)');
ylabel('g(t)');
MW_makeplotlookbetter(14,[],[12.8,4]);
xlim([min(f)-3,max(f)+3]);
ylim([min(g)-3,max(g)+20]);

% plot scatter 2
%hdelayscatter=figure(7); clf; hold on;
subplot(1,3,3); hold on;
%plot(f(n),g,'.','LineWidth',2,'Color',mylinecolors(1,:),'MarkerSize',3^2);
plot(f(nprime-DELTA),g,'.','LineWidth',2,'Color',mylinecolors(4,:),'MarkerSize',3^2);

% cosmetic
%legend('f(n) vs g(n)','f(n) vs g(n+\tau)')
xlabel('f(t+tau)');
ylabel('g(t)');
MW_makeplotlookbetter(14,[],[12.8,4]);
xlim([min(f)-3,max(f)+3]);
ylim([min(g)-3,max(g)+20]);

%%
if exist('SAVEPLEASE','var')
    %figure(hdelay); MW_makeplotlookbetter(14,[],[3.8,4]);
    figure(hdelay); MW_makeplotlookbetter(14,[],[12.8,4]);
    saveas(hdelay,[EXPORTFOLDER 'SVG_technical_delayscatter.svg']);
    saveas(hdelay,[EXPORTFOLDER 'TIF_technical_delayscatter.tif']);
    saveas(hdelay,[EXPORTFOLDER 'FIG_technical_delayscatter.fig']);
end

%% Another interesting measure might be mutual information
% To use this, one would need to estimate the PDF for both measured
% parameters. This could be done easily with the function kde2d.

nrbins=2^8; % 256

data1       = [f(nprime);g]';
data2deltaf = [f(nprime-DELTA);g]';

[bandwidth,density,X,Y]=kde2d(data1,nrbins);
%[bandwidth,density,X,Y]=kde2d(data2deltaf,nrbins);

% This density seems already normalized
dx=X(1,2)-X(1,1);
dy=Y(2,1)-Y(1,1);
%{
sum(density(:).*dx.*dy)
%}

% note: could use X and Y to generate bins for below..

%densityFlipped=flip(density ,1);
%imagesc([X(1,1),X(1,end)],[Y(end,1), Y(1,1)],densityFlipped);
%set(gca,'YDir','normal')

figure; clf;
imagesc([X(1,1),X(1,end)],[Y(1,1), Y(end,1)],density);
set(gca,'YDir','normal')

% using histograms:
figure; clf; 
[n,edges]=histcounts(data1(:,1),nrbins);
centers=(edges(2:end)+edges(1:end-1))/2;
dc=edges(2)-edges(1);
x1=centers;
pdf1=n./(sum(n)*dc);
plot(x1, pdf1);


figure; clf;
[n,edges]=histcounts(data1(:,2),nrbins);
centers=(edges(2:end)+edges(1:end-1))/2;
dc=edges(2)-edges(1);
x2=centers;
pdf2=n./(sum(n)*dc);
plot(x2, pdf2);

%%

% Using the kdepdf
figure; clf; 
densityx=sum(density.*dy);
x1=X(1,:);
pdf1=densityx;
warning('Double check whether pdf matches x now!');
plot(x1, pdf1);

figure; clf; 
densityy=sum(density'.*dx);
x2=Y(:,1);
pdf2=densityy;
warning('Double check whether pdf matches x now!');
plot(x2, pdf2);

%{
figure; clf;
[n,edges]=histcounts(data1(:,2),nrbins);
centers=(edges(2:end)+edges(1:end-1))/2;
dc=edges(2)-edges(1);
x2=centers;
pdf2=n./(sum(n)*dc);
plot(x2, pdf2);
%}
%% p(x)p(y) distribution
figure; clf;

pdf1rep=repmat(pdf1,[size(pdf2,2),1]);
pdf2rep=repmat(pdf2',[1, size(pdf1,2)]);

pdf12=pdf1rep.*pdf2rep;

imagesc(pdf12);
set(gca,'YDir','normal')

%% 

% Now calculate the mutual information
% See: http://mathworld.wolfram.com/MutualInformation.html

% mutual information is defined as 
% I(X; Y) = sum_X sum_Y { P(x,y) log2 P(x,y)./(P(x)P(y)) }
%
% Note that the last term is just the joined probability

% Without the sum, this is basically the ratio between the pdf and the
% multiplied pdfs x and y, weighed by the pdf. This can also be plotted:

mutualInformationMatrix = density.*dx.*dy .* log2(  density.*dx.*dy./(pdf1rep.*pdf2rep.*dx.*dy)  );

myEpsilon = 10e-10;
zeroMask = zeros(size(density));
zeroMask(density.*dx.*dy<myEpsilon)=1;
zeroMask(pdf1rep.*dx<myEpsilon)=1;
zeroMask(pdf2rep.*dy<myEpsilon)=1;

mutualInformationMatrixPrime=mutualInformationMatrix;
mutualInformationMatrixPrime(logical(zeroMask))=0;

%{
% Making entries zeroes also doesn't work that well because some imaginary
numbers remain
mutualInformationMatrixPrime = mutualInformationMatrix;
mutualInformationMatrixPrime(isnan(mutualInformationMatrixPrime)) = 0;
mutualInformationMatrixPrime(isinf(mutualInformationMatrixPrime)) = 0;
%}

%{
% Replace almost zero values with zero values
myEpsilon = 10e-15;
densityEpsilon = density;
densityEpsilon(densityEpsilon<myEpsilon)=0;
pdf1repEpsilon = pdf1rep;
pdf1repEpsilon(pdf1repEpsilon<myEpsilon)=0;
pdf2repEpsilon = pdf2rep;
pdf2repEpsilon(pdf2repEpsilon<myEpsilon)=0;
% And calculate again
mutualInformationMatrixEpsilon = densityEpsilon .* log(  densityEpsilon./(pdf1repEpsilon.*pdf2repEpsilon)  )./log(2);
% NOTE: This doesn't work so well becaue it still produces Inf and NaN
values
%}


figure(10); clf; hold on;

subplot(1,5,1); 
imagesc(density);
set(gca,'YDir','normal')
title('P(x,y)')

subplot(1,5,2); 
imagesc(pdf1rep.*dx);
set(gca,'YDir','normal')
title('P(x)')

subplot(1,5,3); 
imagesc(pdf2rep.*dy);
set(gca,'YDir','normal')
title('P(y)')

subplot(1,5,4); 
imagesc(pdf1rep.*pdf2rep.*dx.*dy);
set(gca,'YDir','normal')
title('P(x)*P(y)')

subplot(1,5,5); 
imagesc( mutualInformationMatrixPrime );
set(gca,'YDir','normal')
title('Terms in sum for I(X;Y)')

toSumOver = mutualInformationMatrixPrime(:);
TheMutualInformation = sum(toSumOver)

%%

% But it's easier using the definition in terms of entropy
% ===

% Note that each entry can be interpreted as a state, so the pdfs should
% be normalized to 1.
pdf1norm = pdf1./sum(pdf1);
pdf2norm = pdf2./sum(pdf2);
densitynorm = density./sum(density(:));

% entropy pdf(x)
H1=-sum(pdf1norm.*log2(pdf1norm))
% entropy pdf(y)
H2=-sum(pdf2norm.*log2(pdf2norm))
% entropy pdf(x,y)
myEpsilon=10-10; % to handle epsilon values
densityEpsilon=densitynorm;
densityEpsilon(densityEpsilon<myEpsilon)=0;
densityTerms=densityEpsilon.*log2(densityEpsilon); % actual calculation of terms
densityTerms(isinf(densityTerms))=0;
densityTerms(isnan(densityTerms))=0;
H12=-sum(sum(densityTerms)) % summation of terms

mutualInformation = H1+H2-H12;

mutualInformationMW = entropyMW(pdf1)+entropyMW(pdf2)-entropyMW(density)

% Refs 
% ===
% https://stackoverflow.com/questions/22074941/shannons-entropy-calculation
% https://nl.mathworks.com/help/images/ref/entropy.html
    % Note that matlab's function "entropy" determines the entropy for a group
    % of observations, i.e. the input is not a pdf function, a pdf function is
    % generated based on the input using imhist.
% http://mathworld.wolfram.com/Entropy.html
% http://mathworld.wolfram.com/MutualInformation.html

% Conversely, we can used the continuous definition of entropy
% ====

% This is:
% h(x) = -Int { f(x) log (f(x)) dx }

% Refs
% ===
% https://en.wikipedia.org/wiki/Differential_entropy

disp(['Also, mutual info is I(X;Y) = H1+H2-H12 = ' sprintf('%0.5f + %0.5f - %0.5f = %0.5f',  H1, H2, H12, mutualInformation)]);

% Note also that log(a/b) = log(a)-log(b).

%% Comparison / sanity check

% Now, if we'd compare p(x)p(y) with p(x)p(y), instead of p(x)p(y) with
% p(x,y), we'll see that mutual information is 0.

% Note that each entry can be interpreted as a state, so the pdfs should
% be normalized to 1.
pdf1norm = pdf1./sum(pdf1);
pdf2norm = pdf2./sum(pdf2);
pdf12mat = pdf1rep .* pdf2rep;
pdf12matNorm = pdf12mat ./ sum(pdf12mat(:));

% entropy pdf(x)
H1=-sum(pdf1norm.*log2(pdf1norm))
% entropy pdf(y)
H2=-sum(pdf2norm.*log2(pdf2norm))
% entropy pdf(x,y)
myEpsilon=10-10; % to handle epsilon values
pdf12matEpsilon=pdf12matNorm;
pdf12matEpsilon(pdf12matEpsilon<myEpsilon)=0;
pdf12matTerms=pdf12matEpsilon.*log2(pdf12matEpsilon); % actual calculation of terms
pdf12matTerms(isinf(pdf12matTerms))=0;
pdf12matTerms(isnan(pdf12matTerms))=0;
H12=-sum(sum(pdf12matTerms)); % summation of terms

mutualInformation = H1+H2-H12;

disp(['For pdf(x)pdf(y), mutual info is I(X;Y) = H1+H2-H12 = ' sprintf('%0.5f + %0.5f - %0.5f = %0.5f',  H1, H2, H12, mutualInformation)]);

