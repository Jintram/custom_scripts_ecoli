
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

if ~ishandle(hratio); hratio=figure(); else figure(hratio); end
clf; hold on;

for i=1:3
    plot([alldata(i).frameData(:).frameTime],fluorValuesAll{i}(1,:)./fluorValuesAll{i}(2,:),'LineWidth',3)
end

title('Ratio between fluor values')
xlabel('time [min]');
ylabel(['mCerulean / mVenus' 10 'or prrsa / pn25']);
MW_makeplotlookbetter(20);