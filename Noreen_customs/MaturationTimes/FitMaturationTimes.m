function [matrate,Fimmat,matrate_confint,Fimmat_confint,mattime,Ftotal]=FitMaturationTimes(lineagenumber, xdata, ydata,  ...
    timeOnsetAB, timeMaxFit, timerangeAllFluoMature,SHOWPLOT)

% input:
% lineagenumber:  any number to identify later on which lineage has been plotted
%           (e.g. 1)
% ydata:   one time traces (which can consist fo data of one branch or the sum of
%           many branches.  ydata should be smoothed and bleaching corrected
%               = total fluorescence data
% xdata:   corresponding time data points
% timeOnsetAB: time at which antibiotics were added (in min), respectively the
%          time after which an exp saturation curve of fluorescence is observed
%          (AB's need some time to diffuse into mebmrane and to cells
% timeMaxFit: max time in fit range (purpose: exluce the very large t when
%          fluorescnce is cst and no more maturation occurs
% timerangeAllFluoMature: [t1 t2]  time window over which the total fluorescence
%               is averaged (to obtain the longterm=all fluorophores are
%               mature plateau value
% SHOWPLOT: 0 or 1: display figures
%
% Maturation process:
% ydata = Fmat + Fimmat * { 1 - exp[-matrate*(xdata-timeonsetAB)]}
% used formula : 
% ydata = Ftotal - Fimmat * exp[-matrate*(xdata-timeonsetAB)]
%       with:  Ftotal = Fimmat + Fmat  (immature and mature total protein amount)
% 
% yforfit = Fimmat *  exp[-matrate*xforfit]
%       with: yforfit = Ftotal - ydata
%             xforfit = xdata - timeonsetAB
% fit parameters are Fimmat (->a in fit function) and matrate (-> b in fit
% function)
%
% fit function: y=a*exp(b*x)
% a-> total amount of immature protein (at start of fit, i.e. supposedly at
% addition of AB's if it was instantaneous)
% b-> maturation rate in [min]

% **********************************************
% *** PREPARE DATA ***



if length(xdata)~=length(ydata)
    error('xdata and ydata must have same length.')
end

% total fluorescence after maturation: Ftotal
idxallmature=find(xdata>timerangeAllFluoMature(1) & xdata<timerangeAllFluoMature(2));
Ftotal = mean(ydata(idxallmature));

% get timewise subsets of xdata and ydata and convert to new fit variables.
% make sure that all vectors are column vectors
if size(xdata,1)==1
    xdata=xdata';
end
if size(ydata,1)==1
    ydata=ydata';
end
idxfittime=find(xdata>timeOnsetAB & xdata<timeMaxFit);
xdatasub=xdata(idxfittime);  % 'sub' -> within fittimerange
ydatasub=ydata(idxfittime);
yforfit=Ftotal-ydatasub;
xforfit=xdatasub-timeOnsetAB;

% **********************************************
% *** PERFORM FIT ***

% (basically copied form "GenerateCode" of Fitting Toolbox)
% check for Nan and Inf values
ok_ = isfinite(xforfit) & isfinite(yforfit);
if ~all( ok_ )
    warning( 'FitMaturationTimes: IgnoringNansAndInfs Ignoring NaNs and Infs in data.' );
end
% suitable start paramters for [a b]
% st_ = [0.94164894964346069 -0.055749474676731262 ];
st_= [ 100000 -0.01];
ft_ = fittype('exp1');   % y = a*exp(bx)  (b is negative!)
% Fit this model using new data
cf_ = fit(xforfit(ok_),yforfit(ok_),ft_,'Startpoint',st_);
% -- likely obsolete start ---  
% Alternatively uncomment the following lines to use coefficients from the
% original fit. You can use this choice to plot the original fit against new
% data.
%    cv_ = { 1.2542353308339842, -0.063780182331307986};
%    cf_ = cfit(ft_,cv_{:});
% --- likely obsolete end ---

% **********************************************
% *** GET FITTED PARAMETERS ***
Fimmat=cf_.a;
matrate=-cf_.b;
mattime=1/matrate;
confinterval=confint(cf_,0.95); % 95% confidence interval (=default)
Fimmat_confint=confinterval(:,1);
matrate_confint=confinterval(:,2);

% **********************************************
% *** SHOW FIGURES ***
if SHOWPLOT

    % fig1: direct fit
    fig1=figure(lineagenumber+10);
    set(fig1,'WindowStyle','docked')
    clf 
    subplot(2,1,1)
    hold on
    plot(xforfit,yforfit,'.b','MarkerSize',15)
    plot(cf_,'r','predfunc',0.95)
    xlabel('xforfit = time - timeOnsetAB')
    ylabel('yforfit=totalFluo-ydata')
    titlestring=['matrate=' num2str(matrate) ...
        '; mattime= ' num2str(mattime) '; 95% conf. interval.' ];
    title(titlestring)

    % fig2: overview
    % create line of fitted function
    xdatafitted=xdatasub(1):1:xdatasub(end);
    ydatafitted=Ftotal - Fimmat * exp(-matrate*(xdatafitted-timeOnsetAB));
   

    %fig2=figure();
    %set(fig2,'WindowStyle','docked')
    %clf 
    subplot(2,1,2)
    hold on
    plot(xdata,ydata,'.b','MarkerSize',15)
    plot(xdatafitted,ydatafitted,'-r')

    currxlim=get(gca,'xlim');
    xlim([timeOnsetAB-40  currxlim(2)+10]);
    % plot time ranges
    currylim=get(gca,'ylim');
    plot([ timeOnsetAB timeOnsetAB],currylim,'k','LineWidth',2)
    plot([ timeMaxFit timeMaxFit],currylim,'k','LineWidth',2)
    plot([ timerangeAllFluoMature(1) timerangeAllFluoMature(1)],currylim,'c','LineWidth',1)
    plot([ timerangeAllFluoMature(2) timerangeAllFluoMature(2)],currylim,'c','LineWidth',1)
    xlabel('time [min]')
    ylabel('total fluo')
    titlestring=['dataset ' num2str(lineagenumber)];
    title(titlestring);
end