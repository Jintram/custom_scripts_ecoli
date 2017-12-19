function [pdf1, pdf2, pdf1pdf2, pdf12, dx, dy] = getpfds12MW(data,nrbins,plottingyesno)
% [pdf1, pdf2, pdf12] = getpfds12MW(data1,nrbins)
% 
% IN
%   data      [x, y], series of observations in stochastic parameters X and
%               Y, N by 2 matrix
%   nrbins    optional nr of bins, default is 2^8 or 256 (needs to be power
%               of 2)
%
% OUT
%   pdf1      distribution p(x)dx
%   pdf2      distribution p(y)dy
%   pdf1pdf2  distribution p(x)p(y)dxdy
%   pdf12     joint distribution p(x,y)dxdy
%   dx        p(x)dx should integrate to 1
%   dy        p(y)dy should integrate to 1
% 
% Note that all output is normalized to 1, such that the output actually 
% equals p(x)dx, p(x,y)dxdy etc.

%% params
if ~exist('nrbins','var')
    nrbins=2^8; % 256
end

if ~exist('plottingyesno','var')
    plottingyesno=0; % 256
end

%%

% use kernel density estimation to get probability distribution
[bandwidth,pdf12,X,Y]=kde2d(data,nrbins);

% calculate the width of the bins
dx=X(1,2)-X(1,1);
dy=Y(2,1)-Y(1,1);

% calculate the marginal distributions
pdf1=sum(pdf12.*dy);
pdf2=sum(pdf12'.*dx);

% calculate p(x)p(y) 
pdf1rep=repmat(pdf1,[size(pdf2,2),1]);
pdf2rep=repmat(pdf2',[1, size(pdf1,2)]);
pdf1pdf2=pdf1rep.*pdf2rep; 

% now normalize all vectors to one
% One way to do this is:
pdf1prime        = pdf1.*dx;
pdf2prime        = pdf2.*dy;
pdf1pdf2prime    = pdf1pdf2.*dx.*dy;
pdf12prime       = pdf12.*dx.*dy;

% an easier one is:
pdf1        = pdf1./sum(pdf1(:));
pdf2        = pdf2./sum(pdf2(:));
pdf1pdf2    = pdf1pdf2./sum(pdf1pdf2(:));
pdf12       = pdf12./sum(pdf12(:));

    
xcorners = [X(1), X(end)];
ycorners = [Y(1), Y(end)];

% 
if plottingyesno
    %% Plotting of the pdfs
    myColors=linspecer(4);
    figure(1); clf; hold on;
    
    subplot(1,4,1); hold on;
    plot(X(1,:),pdf1,'LineWidth',2,'Color',myColors(1,:));
    plot(X(1,:),pdf1prime,'--','LineWidth',2,'Color',myColors(2,:));
    MW_makeplotlookbetter(16);
    xlabel('x');
    ylabel('p(x)');
    xlim([min(X(1,:)),max(X(1,:))]);
    
    subplot(1,4,2); hold on;
    plot(Y(:,1),pdf2,'LineWidth',2,'Color',myColors(3,:));
    plot(Y(:,1),pdf2prime,'--','LineWidth',2,'Color',myColors(4,:));
    MW_makeplotlookbetter(16);
    xlabel('y');
    ylabel('p(y)');
    xlim([min(Y(:,1)),max(Y(:,1))]);
    
    subplot(1,4,3); hold on;
    %axis off
    imagesc(xcorners,ycorners,pdf1pdf2)
    set(gca,'YDir','normal')
    xlim(xcorners)
    ylim(ycorners)
    %xlim([1,size(pdf1pdf2,2)])
    %ylim([1,size(pdf1pdf2,1)])
    xlabel('x');
    ylabel('y');    
    MW_makeplotlookbetter(16);
    
    subplot(1,4,4); hold on;
    %axis off
    imagesc(xcorners,ycorners,pdf12)
    xlim(xcorners)
    ylim(ycorners)
    %xlim([1,size(pdf12,2)])
    %ylim([1,size(pdf12,1)])
    xlabel('x');
    ylabel('y');
    MW_makeplotlookbetter(16);
    
    set(gca,'YDir','normal')
end








