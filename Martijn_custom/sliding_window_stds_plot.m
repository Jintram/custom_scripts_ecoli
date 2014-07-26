function [lgh1, lgh2, lgh3, lgh4] = sliding_window_stds_plot(myPhosphoData,myPhosphoAuxiliary,datax,datay,fignumstart,mylinecolor,showplots,divbymean,ID)
    % Set nrybins and nrtbins in this function!
    %
    % lgh1, lgh2, lgh3 are line handles, which can be used to make a
    % legend.

    % obtain statistics per time window
    nrybins=20; %30
    nrtbins=1; %10
    [t_center,means,stds,yvalue_n,lgh1] = sliding_window_stds(nrybins,datax,datay,nrtbins,showplots,fignumstart,mylinecolor);

    % Restcale stds by means
    if divbymean
        stds = stds./means;
    end
    
    % make fits
    means_fitCoef1 = polyfit(t_center,means,1);
    means_fitted = means_fitCoef1(1)*t_center + means_fitCoef1(2);    
    
    std_fitCoef1 = polyfit(t_center,stds,1);
    std_fitted = std_fitCoef1(1)*t_center + std_fitCoef1(2);
    
    x_values=[0:.01:2];
    stdmean_fitCoef1 = polyfit(means,stds,1);
    stdmean_fitted = stdmean_fitCoef1(1)*x_values + stdmean_fitCoef1(2);
       
    % Plot means per timewindow
    % ===
    
    figure(fignumstart+1); hold on;
    lgh2 = plot(t_center,means,'-o','LineWidth',3,'color',mylinecolor);
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
    lgh3 = plot(t_center,stds,'-o','LineWidth',3,'color',mylinecolor);
    plot(t_center,std_fitted,'-','LineWidth',3,'color','k');
    
    oldylim = ylim; newylim = [0, max(stds)*1.1];
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
    
    % Plot standard deviations vs. mu
    % ===
    
    figure(fignumstart+3); 
    hold on;
    lgh4 = plot(means,stds,myPhosphoAuxiliary.markers.(ID),'LineWidth',3,'color',mylinecolor);    
    %plot(x_values,stdmean_fitted,'-','LineWidth',3,'color','k');
    
    axis([0,2,0,2])
    
    set(gca,'FontSize',20);
    xlabel('Growth speed (dlb/hr)');    
    if divbymean
        ylabel(['Standard deviation / mean']);    
        title('Growth speed vs noise');
    else
        ylabel(['Standard deviation']);    
        title('Growth speed vs std. dev.');
    end;
    

end