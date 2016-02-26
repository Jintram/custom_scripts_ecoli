

%%

PATHNAME = 'F:\A_Tans1_step1_incoming_not_backed_up\2016_02-01_flowcellMWtest\pos1\';
FILEBASE = 'pos1-g-';

%%

meanFluorValues = [];
times = [];
for i = 5:1325

    myfilepath = [PATHNAME FILEBASE sprintf('%03d', i) '.tif'];
    
    theimg = imread(myfilepath);
    
    myinfo=imfinfo(myfilepath);
    file_date=myinfo.FileModDate;
    
        
    times(end+1) = datenum(file_date);    
    meanFluorValues(end+1) = mean(theimg(:));
    
end

%% plot

SWITCHTIME = 8;

myYLim = [0,max(meanFluorValues)*1.1];

t0 = times(1);
timesrelt0inmins = (times-t0)*24*60;
mycolors = linspecer(10);

figure(1); clf; hold on;
plot(timesrelt0inmins,meanFluorValues,'Color',mycolors(1,:),'LineWidth',3);
xlabel('Time [mins]');
ylabel('Fluor signal [a.u.]')
title(['Fluoresceine flown into normal gel pad' 10 'script20160201_flowcellMWtest.m'],...
        'Interpreter','None');
MW_makeplotlookbetter(20);

% plot switch
plot([SWITCHTIME,SWITCHTIME],myYLim,':k','LineWidth',2)


%% theoretical prediction
%See Mathematica sheet 2016_02_02_diffusionForPad

% diffusion constant
alpha=1/100;
% time
t = [1:600];
% system
rprime = 2.5*10^-3;
tstar = 600;
% geometrical/time constraints
r=0;
tau = 0;

% function
GR00=...
    1./(4*pi*alpha*(t-tau)).*...
    exp((-(r^2+rprime^2))./(4*alpha*(t-tau)));%.*...
    %besselj(0,(r*rprime)./(2*alpha*(t-tau)));
    % note that besselj falls out since r*rprime is zero, which makes it a
    % simple exponential function.
GR00integral=...
    cumsum(GR00);

% plot
figure(2); clf; hold on;
plot(t,GR00integral,'Color',mycolors(2,:),'LineWidth',3);

xlabel('Time [mins]');
%ylabel('Fluor signal [a.u.]')
%title(['Fluoresceine flown into normal gel pad' 10 'script20160201_flowcellMWtest.m'],...
%        'Interpreter','None');
MW_makeplotlookbetter(20);


%{
GR00[r_, t_, rprime_, tau_] := 
  1/(4*pi*tau (t - tau))*
    Exp[(-(r^2 + rprime^2))/(4*tau*(t - tau))]*
    BesselJ[0, (r*rprime)/(2*tau (t - tau))] /. {tau -> 
     1/2};
GR00integral[r_, tstar_, rprime_, tau_] := 
  NIntegrate[GR00[r, t, rprime, tau]
   , {t, 0, tstar}];
myradius = 2.5*10^-3;
Plot[GR00[0, t, myradius, 0], {t, 0, 600}]
Plot[GR00integral[0, t, myradius, 0], {t, 0, 600}]
%}





