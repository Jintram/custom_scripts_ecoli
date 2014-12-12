% xx and yy are random distribution with either no or a given correlation
% coefficient.
% plots the scatter plot and the binned yy-values. The latter can be used
% to determine the slope=transmission coefficient.
% Question to answer: Is the slope not actually the same as the correlation
% coefficient? see nature paper SOM
%
% Use vectors of dimension nx1, not 1xn

counter=1;
    for i=min(xx):0.3:max(xx)  
         bin= xx>i &  xx<= i+0.3;       
         xbin(counter,1)  = mean(xx(bin)); 
         ybin(counter,1)    = mean(yy(bin));
         counter          = counter+1;
    end

pol=polyfit(xx,yy,1);
myoffset=pol(2);
myslope=pol(1);
    
clf
plot(xx,yy,'r.'); hold on
plot(xbin,ybin,'b.','MarkerSize',15); 
xvec=min(xx):0.01:max(xx);
plot(xvec,xvec*myslope+myoffset,'k')
title([ 'Corrcoeff=' num2str(corr(xx,yy)) ',   Slope (Transmission)= ' num2str(myslope)]);


%% same as above but with zz instead of yy.
% usefull if e.g. yy and xx are uncorrelated and then zz has a
% preset corrcoeff rho with xx:
% zz=rho*yy + sqrt(1-rho^2)*xx

counter=1;
    for i=min(xx):0.3:max(xx)  
         bin= xx>i &  xx<= i+0.3;       
         xbin(counter,1)  = mean(xx(bin)); 
         zbin(counter,1)    = mean(zz(bin));
         counter          = counter+1;
    end

pol=polyfit(xx,zz,1);
myoffset=pol(2);
myslope=pol(1);
    
clf
plot(xx,zz,'r.'); hold on
plot(xbin,zbin,'b.','MarkerSize',15); 
xvec=min(xx):0.01:max(xx);
plot(xvec,xvec*myslope+myoffset,'k')
title([ 'Corrcoeff=' num2str(corr(xx,zz)) ',   Slope (Transmission)= ' num2str(myslope)]);


%%
% the same but with different variable names
clear xbin ybin counter myoffset myslope myxvector myyvector

% ** adjust **
myxvector=conc;%norm;
myyvector=mu;%norm;
myxlabel='conc E';
myylabel='mu';
increment=0.1;
% ***

counter=1;
    for i=min(myxvector):increment:max(myxvector)  
         bin= myxvector>i &  myxvector<= i+increment;       
         xbin(counter,1)  = mean(myxvector(bin)); 
         ybin(counter,1)    = mean(myyvector(bin));
         counter          = counter+1;
    end

pol=polyfit(myxvector,myyvector,1);
%pol=polyfit(xbin,ybin,1);
myoffset=pol(2);
myslope=pol(1);
    
clf
plot(myxvector,myyvector,'r.'); hold on
plot(xbin,ybin,'b.','MarkerSize',15); 
xvec=min(myxvector):0.01:max(myxvector);
plot(xvec,xvec*myslope+myoffset,'k')
title([ 'Corrcoeff=' num2str(corr(myxvector,myyvector)) ',   Slope (Transmission)= ' num2str(myslope)]);
xlabel(myxlabel);
ylabel(myylabel);
