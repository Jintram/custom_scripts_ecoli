function sliding_window_stds_plot(datax,datay,fignumstart,mylinecolor,showplots,divbymean)
    % Set nrybins and nrtbins in this function!

    % obtain statistics per time window
    nrybins=50;
    nrtbins=10;
    [t_center,means,stds,yvalue_n] = sliding_window_stds(nrybins,datax,datay,nrtbins,showplots,fignumstart,mylinecolor);

    % Restcale stds by means
    if divbymean
        stds_rescaled = stds./means;
    end
    
    % make fits
    means_fitCoef1 = polyfit(t_center,means,1);
    means_fitted = means_fitCoef1(1)*t_center + means_fitCoef1(2);    
    
    std_fitCoef1 = polyfit(t_center,stds_rescaled,1);
    std_fitted = std_fitCoef1(1)*t_center + std_fitCoef1(2);
       
    % Plot means per timewindow
    % ===
    
    figure(fignumstart+1); hold on;
    plot(t_center,means,'-o','LineWidth',3,'color',mylinecolor);
    plot(t_center,means_fitted,'-','LineWidth',3,'color','k');
    
    oldylim = ylim; newylim = [0, max(means)*1.1];
    if oldylim(2) < newylim(2) % set new ylim if adjustment needed
        ylim(newylim);
    end
    set(gca,'FontSize',20);
    title('Mean \mu (dbl/hr)');
    xlabel('Points in time');
    ylabel('Mean \mu (dbl/hr)');
    
    % Plot standard deviations per timewindow
    % ===
    
    figure(fignumstart+2); hold on;
    plot(t_center,stds_rescaled,'-o','LineWidth',3,'color',mylinecolor);
    plot(t_center,std_fitted,'-','LineWidth',3,'color','k');
    
    oldylim = ylim; newylim = [0, max(stds_rescaled)*1.1];
    if oldylim(2) < newylim(2) % set new ylim if adjustment needed
        ylim(newylim);
    end
    set(gca,'FontSize',20);
    xlabel('Points in time');
    % if div. mean, show in label
    if divbymean
        title('Noise');
        ylabel(['Standard deviation / mean']);    
    else 
        title('Standard deviations');
        ylabel(['Standard deviation']);        
    end;
    

end