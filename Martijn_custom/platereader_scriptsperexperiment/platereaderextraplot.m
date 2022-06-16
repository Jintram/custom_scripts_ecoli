function [allLinehandles, legendLinesHandles] = platereaderextraplot(USERSETTINGS,mydata,myhandle,legendLinesHandles,mycolor)
    
    figure(myhandle);
    
    % y=zero line
    plot(USERSETTINGS.myXlim, [0, 0], '-','LineWidth',2,'Color','k')
    
    allLinehandles=[];
    for i = 1:numel(mydata)
        l=plot(   ones(1,numel(mydata(i).manualMuValues)) * USERSETTINGS.myConcentrationValues(i), ...
                    mydata(i).manualMuValues,'o')
        set(l, 'MarkerSize', 10,'LineWidth',3,'Color',mycolor);   
        allLinehandles(end+1)=l;
    end
    legendLinesHandles(end+1) = l;
    % value at zero
    plot(USERSETTINGS.myXlim,[1,1]*mean(mydata(end).manualMuValues),'--','LineWidth',3,'Color',mycolor);

    % Some settings for the plot
    set(gca,'xscale','log');

    MW_makeplotlookbetter(20);

    ylim([-.1,1]);
    xlim(USERSETTINGS.myXlim);

    xlabel('cAMP concentration [M]');
    ylabel('fitted growth rate [dbl/hr]');
end