




%% 

% This script has been moved to 
% script20160429_filamentRecoveryDivisionRatioss





















%% 
% Plotting division ratios
%
% This script requires the "schnitzcells" struct to exist and 

%% 
FIGURE=3;
LENGTHFIELD = 'areaPixels';
%LENGTHFIELD = 'length_fitNew';
%LENGTHFIELD = 'length_skeleton'

%%

myLengthNewborns = []; myLengthParents = [];
for i = 1:numel(schnitzcells)
    
    LengthNewborn = schnitzcells(i).(LENGTHFIELD)(1);
    
    parentSchnitz = schnitzcells(i).P;    
    
    if parentSchnitz ~=0
    
        LengthParent = schnitzcells(parentSchnitz).(LENGTHFIELD)(end);
        
        myLengthNewborns(end+1) =   LengthNewborn;
        myLengthParents(end+1) =    LengthParent;
        
    end
    
end

%%
% user given parameters
LONGESTNORMALDIVSIZEPARENT = 2000;
LEFTX =  0;
RIGHTX = 30;

% calculated parameters
if strcmp(LENGTHFIELD,'length_skeleton')
    rightx = RIGHTX;
else
    rightx  = max(myLengthParents)*1.1;
end

% create figure
figure(FIGURE); clf; hold on;

% plot helping lines at 1/2n
% for i=1:5
%     plot([rightx, LEFTX], [.5/i .5/i],':','Color',[.5 .5 .5],'LineWidth',2)
%     plot([rightx, LEFTX], 1-[.5/i .5/i],':','Color',[.5 .5 .5],'LineWidth',2)
% end
N=5;
for i=1:N
    for j = 1:(i*2-1)
        plot([0, rightx], [(j)/(2*i) (j)/(2*i)],'-','Color',[.5 .5 .5],'LineWidth',N-i+1)
    end
end

% plot ratios
Ratios = myLengthNewborns./myLengthParents;
plot(myLengthParents,Ratios,'xb','LineWidth',2);
plot(myLengthParents,1-Ratios,'xr','LineWidth',2);

% target length line
x = (LONGESTNORMALDIVSIZEPARENT):rightx;
y = (LONGESTNORMALDIVSIZEPARENT/2)./x;
plot(x,y,'-','Color',[.5 .5 .5],'LineWidth',2);



% cosmetics
ylim([0,1]);
xlim([LEFTX,rightx]);

xlabel(['Cell size by ' LENGTHFIELD],'Interpreter','None');
ylabel('L_{child}/L_{parent}');

MW_makeplotlookbetter(15)













