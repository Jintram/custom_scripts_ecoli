
%% 
% Plotting division ratios
%

%% The sulA datasets

datasetsPaths = ...
    { ...
    'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos1crop\data\pos1crop-Schnitz.mat',...
    'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos4crop\data\pos4crop-Schnitz.mat',...
    'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-08_sulA_recovery_200uM_IPTG\pos7crop\data\pos7crop-Schnitz.mat',...
    }



for datasetIdx = 1:3 
    
    load(datasetsPaths{datasetIdx});

    %% 
    FIGURE=datasetIdx;
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

    % calculated parameters
    if strcmp(LENGTHFIELD,'length_skeleton')
        rightx = 30;
    elseif strcmp(LENGTHFIELD,'areaPixels')
        rightx = 20000;
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

end











