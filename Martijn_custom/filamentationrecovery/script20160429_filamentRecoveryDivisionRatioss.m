
%% 
% Plotting division ratios
%

HISTNRBINS=20;

FIGURENUMBERS=[1 1 1];

LENGTHFIELD = 'areaPixels';
%LENGTHFIELD = 'length_fitNew';
%LENGTHFIELD = 'length_skeleton';

%% The sulA datasets

datasetsPaths = ...
    { ...
    'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos1crop\data\pos1crop-Schnitz.mat',...
    'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos4crop\data\pos4crop-Schnitz.mat',...
    'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos7crop\data\pos7crop-Schnitz.mat',...
    }

%% 
for i = FIGURENUMBERS
    figure(i); clf; hold on;
end

%%

for datasetIdx = 1:numel(datasetsPaths)
    
    load(datasetsPaths{datasetIdx});

    %% 
    figureIndex=FIGURENUMBERS(datasetIdx);

    %%

    
    % Finding parent with each daughter
    %{
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
    %}
    
    % Finding 2 daughters with each parent
    myLengthNewborns{datasetIdx} = []; myLengthParents{datasetIdx} = []; myLengthSumNewborns{datasetIdx}=[];
    for i = 1:numel(schnitzcells)

        LengthParent = schnitzcells(i).(LENGTHFIELD)(end);
        
        daughterSchnitz1 = schnitzcells(i).D;
        daughterSchnitz2 = schnitzcells(i).E;

        if ~any([daughterSchnitz1,daughterSchnitz2]==0)

            lengthDaughterSchnitz1 = schnitzcells(daughterSchnitz1).(LENGTHFIELD)(1);
            lengthDaughterSchnitz2 = schnitzcells(daughterSchnitz2).(LENGTHFIELD)(1);            
            
            % daughter 1
            myLengthNewborns{datasetIdx}(end+1) =     lengthDaughterSchnitz1;
            myLengthParents{datasetIdx}(end+1) =      LengthParent;
            myLengthSumNewborns{datasetIdx}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2;
            % daughter 2
            myLengthNewborns{datasetIdx}(end+1) =     lengthDaughterSchnitz2;
            myLengthParents{datasetIdx}(end+1) =      LengthParent;
            myLengthSumNewborns{datasetIdx}(end+1) =  lengthDaughterSchnitz1+lengthDaughterSchnitz2;

        end

    end

end

%% Now plot this
Ratios={};
for datasetIdx = 1:numel(datasetsPaths)
    
    %%
    % user given parameters
    LONGESTNORMALDIVSIZEPARENT = 2000;
    LEFTX =  0;   

    % calculated parameters
    if strcmp(LENGTHFIELD,'length_skeleton')
        rightx = 30;
    elseif strcmp(LENGTHFIELD,'areaPixels')
        rightx = 20000;
    else
        rightx  = max(myLengthParents)*1.1;
    end

    % create figure
    figure(figureIndex); 

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

    % target length line
    x = (LONGESTNORMALDIVSIZEPARENT):rightx;
    y = (LONGESTNORMALDIVSIZEPARENT/2)./x;
    plot(x,y,'-','Color',[.5 .5 .5],'LineWidth',2);    
    
    % plot ratios
    Ratios{datasetIdx} = myLengthNewborns{datasetIdx}./myLengthSumNewborns{datasetIdx};
    plot(myLengthSumNewborns{datasetIdx},Ratios{datasetIdx},'xb','LineWidth',2);
    plot(myLengthSumNewborns{datasetIdx},1-Ratios{datasetIdx},'ob','LineWidth',2);

    % cosmetics
    ylim([0,1]);
    xlim([LEFTX,rightx]);

    xlabel(['Cell size by ' LENGTHFIELD],'Interpreter','None');
    ylabel('L_{child}/L_{parent}');

    MW_makeplotlookbetter(15)

end

%% Now calculate overall histogram
h=figure(FIGURENUMBERS(end)+1); clf; hold on;

% helping lines
N=5; highestcount=max(count);
for i=1:N
    for j = 1:(i*2-1)
        plot([(j)/(2*i) (j)/(2*i)],[0, highestcount],'-','Color',[.5 .5 .5],'LineWidth',N-i+1)
    end
end

% hist
[count,x] = hist([Ratios{:}],HISTNRBINS);

if exist('histSkel','var') && exist('histArea','var')
    % plot histograms from multiple fields
    l1=plot(histArea(2,:),histArea(1,:),'-','LineWidth',3);    
    l2=plot(histSkel(2,:),histSkel(1,:),'-','LineWidth',3);
    legend([l1,l2],{'Skeleton','Area'});
else
    % plotting of one histogram
    plot(x,count,'-','LineWidth',3);
end

% additional ticks
plot([0:.1:1],zeros(1,11),'+','Color','k');%,'MarkerFaceColor','k');

% cosmetics
%ax=gca; ax.XTick = [0:.1:1];
ylabel('Count');
xlabel('L_d/L_p');

set(gca,'XTick',[0:.1:1])

ylim([0 highestcount]);

MW_makeplotlookbetter(15);

% Save the histogram
if strcmp(LENGTHFIELD, 'areaPixels')
    histArea=[count;x]
elseif strcmp(LENGTHFIELD, 'length_skeleton')
    histSkel=[count;x]
end

%% Sanity check
% load(datasetsPaths{2})
figure(FIGURENUMBERS(end)+2); clf; hold on

% x=y line
plot([1,10^6],[1,10^6],'-k')

% data
for i=1:numel(datasetsPaths)
    plot(myLengthSumNewborns{datasetIdx}, myLengthParents{datasetIdx},'x')
end

% cosmetics
axis equal;
xlim([0 max([myLengthSumNewborns{:}])]);
ylim([0 max([myLengthParents{:}])]);



%%
figure(FIGURENUMBERS(end)+3); clf; hold on;

plot(histArea(2,:),histArea(1,:),'-r');
plot(histSkel(2,:),histSkel(1,:),'-b');







