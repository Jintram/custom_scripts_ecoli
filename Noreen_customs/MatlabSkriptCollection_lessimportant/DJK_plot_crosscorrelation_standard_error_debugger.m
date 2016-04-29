function DJK_plot_crosscorrelation_standard_error_store(p,branch_groups, fieldX, fieldY,varargin);
% Daans function but saves the generated plots
% construction of 'saveDir'very akward! (change ...)

% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bias = 0;      % 0 will adjust for less data at larger delay times
weighing = 2;  % 2 performs 3/4 weighing
extraNorm = 0; % 0 performs no extra normalization    %blubb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'DJK_plot_crosscorrelation_standard_error_store';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
  errorMessage = sprintf('%s\n%s',['Error width input arguments of ' functionName],['Try "help ' functionName '".']);
  error(errorMessage);
end

numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf('%s\n%s',['This input argument should be a String: ' num2str(varargin{i})],['Try "help ' functionName '".']);
      error(errorMessage);
    end
    fieldName = DJK_schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
  end
end
%--------------------------------------------------------------------------

% create SaveDir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~existfield(p,'selectionName ')
  p.selectionName = ''
end
% Just in case it has not been set yet
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'schnitzcells' filesep];
end
if length(p.selectionName) > 0
  p.DJK_saveDir = [p.DJK_saveDir p.selectionName filesep];
end  
% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveDir=p.DJK_saveDir;

if ~existfield(p,'plot_lineage')
    p.plot_lineage=[];
end
plot_lineage=p.plot_lineage;


% CALC COMPOSITE CORR FOR EACH GROUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = struct; p.movieName = 'nonsense';
for i = 1:length(branch_groups)
  [trash, branch_groups(i).composite_corr] = DJK_getCrossCor(p, branch_groups(i).branches, fieldX, fieldY, 'bias', bias, 'weighing', weighing, 'extraNorm', extraNorm);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% crosscorr for each single lineage (makes more sense if just one branch is
% regarded)
for i = 1:length(branch_groups)
    for j=1:length(branch_groups(i).branches)
        [trash, branch_groups(i).composite_corr_single(j)] = DJK_getCrossCor(p, branch_groups(i).branches(j), fieldX, fieldY, 'bias', bias, 'weighing', weighing, 'extraNorm', extraNorm);
    end
end



% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Figure Name
figureName1 = ['crosscorr_' fieldX ' _ ' fieldY 'singletraces'];
figureName2 = ['crosscorr_' fieldX ' _ ' fieldY 'errors'];
figureName3 = ['crosscorr_' fieldX ' _ ' fieldY 'errors_norm'];


% actual plotting 
for i = 1:length(branch_groups)
  composite_corr(i,:) = branch_groups(i).composite_corr.Y;
end

figure(5); 
clf
for i = 1:length(branch_groups)
    %c=0.2+0.8*i/length(branch_groups);
    co=[1 0 0; 0 0 1; 0 1 0; 0.8 0.5 0; 1 0 1];
    %plot(branch_groups(i).composite_corr.X/60, branch_groups(i).composite_corr.Y, '-', 'LineWidth', 2, 'Color', [c 0.5*c]); hold on;
    plot(branch_groups(i).composite_corr.X/60, branch_groups(i).composite_corr.Y, '-', 'LineWidth', 2, 'Color', co(i,:)); hold on
end
plot(branch_groups(i).composite_corr.X/60, mean(composite_corr), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscorr');
title(['single traces (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot([-100 100],[0 0],'-k'); plot([0 0 ],[-100 100],'-k');  % plot axis
legend('1','2','3','4')
set(gca,'xlim',x_lim);set(gca,'ylim',y_lim);
saveas(gcf,[saveDir figureName1 '.png'], 'png');

%figure;
%errorbar(branch_groups(i).composite_corr.X/60,mean(composite_corr),std(composite_corr), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
%xlabel('time [h]');
%ylabel('crosscorr');
%title(['errorbars (' fieldX ', ' fieldY ')'],'interpreter','none');
%hold on
%x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
%plot([-100 100],[0 0],'-k'); plot([0 0 ],[-100 100],'-k');  % plot axis
%set(gca,'xlim',x_lim);set(gca,'ylim',y_lim);
%saveas(gcf,[saveDir figureName2 '.png'], 'png');

%figure;
%errorbar(branch_groups(i).composite_corr.X/60,mean(composite_corr),std(composite_corr)/sqrt(length(branch_groups)), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
%xlabel('time [h]');
%ylabel('crosscorr');
%title(['errorbars normalized (' fieldX ', ' fieldY ')'],'interpreter','none');
%hold on
%x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
%plot([-100 100],[0 0],'-k'); plot([0 0 ],[-100 100],'-k');  % plot axis
%set(gca,'xlim',x_lim);set(gca,'ylim',y_lim);
%%saveas(gcf,[saveDir figureName3 '.png'], 'png');

plot_lineage=[10 30 50 ];
figure(7); 
clf
for i = 4 %1:length(branch_groups) blubb
    if isempty(plot_lineage)
       for j=1:length(branch_groups(i).branches)
         c=0.3+0.7*j/length(branch_groups(i).branches);
         plot(branch_groups(i).composite_corr_single(j).X/60, branch_groups(i).composite_corr_single(j).Y, '-', 'LineWidth', 2, 'Color', [c 0.5 0.5]); hold on;
       end
    else
        runner=0;
       for j=plot_lineage
            runner=runner+1;
            c=0.3+0.7*runner/length(plot_lineage);
            plot(branch_groups(i).composite_corr_single(j).X/60, branch_groups(i).composite_corr_single(j).Y, '-', 'LineWidth', 2, 'Color', [c 0.5 0.5]); hold on;
       end
    end
             
end
plot(branch_groups(i).composite_corr.X/60, branch_groups(i).composite_corr.Y, '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscorr');
title(['single traces (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot([-100 100],[0 0],'-k'); plot([0 0 ],[-100 100],'-k');  % plot axis
str1='0';str2='0';str3='0';str4='0';
if isempty(plot_lineage)
    legend('1','2','3','4');
else
    if length(plot_lineage<4) plot_lineage=[plot_lineage,0,0,0]; end
    str1=num2str(plot_lineage(1)); str2=num2str(plot_lineage(2)); str3=num2str(plot_lineage(3)); str4=num2str(plot_lineage(4));
    legend(str1,str2,str3,str4)
end
set(gca,'xlim',x_lim);set(gca,'ylim',y_lim);
%saveas(gcf,[saveDir figureName1 '.png'], 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%