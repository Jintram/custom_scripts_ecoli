
%% Plot per cell
figure;

names={}; alldata = [];
% loop over datasets
for ii = 1:numel(outputfluorvalues)
    
    hold on;
    
    % loop over values each image
    for jj = 1:numel(outputfluorvalues(ii).dirmultiplePercentile09Fluor)
        data=outputfluorvalues(ii).dirmultiplePercentile09Fluor{jj};
        % plot them separately
        plot(ones(size(data))*ii, data,'x','LineWidth',2);                                
    end
        
    % store the names
    names{end+1}=outputfluorvalues(ii).ValuesNames(end-1);
    
    % 
    alldata=[alldata data];
    
end

xlim([0,8]);
ylim([0 max(alldata)*1.2]);

xlabel('Sample')
ylabel(['Fluorescence values of cells' 10 '(Percentile09) [a.u.]'])

set(gca, 'XTick',[1:numel(names)],'XTickLabel',names);

MW_makeplotlookbetter(15);


%% Plot means

figure;

%loop over datasets
for ii = 1:numel(outputfluorvalues)
    
    hold on;
    
    % plot means
    plot(ii, outputfluorvalues(ii).ValuesMeanSignalNoise,'o','LineWidth',4);
        
end

xlim([0,8]);
ylim([0, max([outputfluorvalues(:).ValuesMeanSignalNoise])*1.2]);

xlabel('Sample')
ylabel(['Mean fluorescence cell areas [a.u.]'])

set(gca, 'XTick',[1:numel(names)],'XTickLabel',names);

MW_makeplotlookbetter(15);
