function [branchData, time_lincomb, yfield_lincomb] = NW_plotLineagesTimeTrace (p, schnitzcells,timefield, yfield,extint,extranormalize,nrbranches)
    
% plots the single branches over time. written for maturation time fitting.
% If a parent cell has progeny, the yfield values of the progeny are either
% added (extensive, e.g. G5sum) or averaged (intensive, e.g. G6mean)
%lineage data is obtained via DJK_getBranches 
% **so far only 'added' possible ***
%
% no fool proof function!
%
% REQUIRED ARGUMENTS
% p
% schnitzcells
% timefield : e.g. 'G_time'
% yfield :    e.g. 'G5_sum'
% extint :    'extensive' or 'intensive' SO FAR ONLY EXTENSIVE IMPLEMENTED!
% extranormalize:  =1 all brach-traces are normalized to=1
%                   =0 no extra norm
% nrbranches : how many branches (individual time traces). minimum: initial
%               #of schnitzes (done automatically). maximum: final # of schnitzes
%
%
% OUTPUT (TO COME)
% each row of yfield (e.g. (1,:) ) correponds to one timetrace

if strcmp(extint,'extensive')==0 & strcmp(extint,'extensive')==0
    error('extint must be  ''extensive'' or ''intensive''. ');
end

branchData = DJK_getBranches(p,schnitzcells,'dataFields',{timefield,yfield});
trimmed_branches = DJK_trim_branch_data(branchData,nrbranches);
branchData=trimmed_branches;


% get initial schnitz numbers
initSchnNrs=[];
for i=1:length(branchData)
    initSchnNrs=[initSchnNrs, branchData(i).schnitzNrs(1)];
end

initSchnitzNrs_unique=unique(initSchnNrs); 

% time data points and yfield data points, averaged such that one vector
% corresponds to all lineages descending from one startschnitz
time_lincomb=branchData(1).(timefield);  % 1 x #datapts . vector
yfield_lincomb=zeros(length(initSchnitzNrs_unique),length(time_lincomb));   % #uniqueSchnitzes x # datapts . matrix
% "each row is one combined lineage"

for i=1:length(initSchnitzNrs_unique)
    % combine lineages which have the same parent schnitz
    idx=find(initSchnNrs==initSchnitzNrs_unique(i));
    subbranches=branchData(idx);
    
    % ********** SO FAR ONLY EXTENSIVE IMPLEMENTED!!! ************
    for k=1:length(subbranches)
        yfield_lincomb(i,:)=yfield_lincomb(i,:)+(subbranches(k).(yfield))./(subbranches(k).count);
        % *** here: implement counting for intensive:..
        % +subbranches(k).(yfield) / subbranches(k).(count) (or sth like
        % that)
    end
    
    % ---- Normalize each branch to its average --- (de)activate?
    if extranormalize==1
        yfield_lincomb(i,:)=yfield_lincomb(i,:)./mean(yfield_lincomb(i,:));
    end
    % -----
end

% plot
figure(1)
set(gcf,'WindowStyle','docked')
clf
hold on
for i=1:size(yfield_lincomb,1)
    plot(time_lincomb,yfield_lincomb(i,:),'.-','Color',[0 0.7 i/size(yfield_lincomb,1)])
end
plot(time_lincomb,mean(yfield_lincomb),'.-k','MarkerSize',10,'LineWidth',2)
xlabel([timefield '(min)'],'Interpreter','None')
ylabel(yfield,'Interpreter','None')


figure(2)
set(gcf,'WindowStyle','docked')
clf
hold on
errorbar(time_lincomb,sum(yfield_lincomb),std(yfield_lincomb))
xlabel([timefield '(min)'],'Interpreter','None')
ylabel(['sum of ' yfield],'Interpreter','None')


    





