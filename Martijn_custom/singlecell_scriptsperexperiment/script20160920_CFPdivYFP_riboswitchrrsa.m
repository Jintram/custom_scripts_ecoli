
% executed first:
%{
MYDIR='G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-09-20_asc976normal\outputSummary\';
FLUORIDXTOPLOT = [1 2];
 
CUSTOMCOLORS=linspecer(3);

SWITCHTIME=119;
FITWINDOWS = [0,SWITCHTIME; SWITCHTIME,800];

MW_summaryplotspreliminaryanalysis
MW_summaryplotspreliminaryanalysisforswitch
%}

if ~exist('hratio','var'), hratio=figure(); elseif ~ishandle(hratio); hratio=figure(); figure(hratio); end 
clf; hold on;

datapile=[];
for i=1:3
    xdata=[alldata(i).frameData(:).frameTime];
    ydata=[fluorValuesAll{i}(1,:)./fluorValuesAll{i}(2,:)];
    plot(xdata(~isnan(ydata)),...
         ydata(~isnan(ydata)) ,'-x','LineWidth',3)
    datapile=[datapile ydata(~isnan(ydata))];
end

theYlim=[0, max(datapile)*1.2];
ylim(theYlim)

%plot switchtime
plot([SWITCHTIME,SWITCHTIME],theYlim,':k','LineWidth',3);

title('Ratio between fluor values')
xlabel('time [min]');
ylabel(['mCerulean / mVenus' 10 'or prrsa / pn25']);
MW_makeplotlookbetter(20);

outputFilepath = [MYDIR 'summaryplot_ceruleanDivVenus'];
saveas(hratio, [outputFilepath '.tif']);
saveas(hratio, [outputFilepath '.fig']);

