
% Plot total colony length vs. time (for multiple colonies)
% ===
%
% This script uses DJK_analyzeMu to obtain timetrace of total colony length
% and then plots this. It uses a struct with multiple schnitzcell
% structures as input.
%
% MW 2014/12
%
%




% Loop over strains contained by myPhosphoData
strains = fieldnames(myPhosphoData)
for str_idx = 1:numel(strains)

    current_strain = char(strains(str_idx));
    repetitions = fieldnames(myPhosphoData.(current_strain))
    
    % Loop over repetitions done for that strain
    for rep_idx = 1:numel(repetitions)
    
        current_repetition = char(repetitions(rep_idx));        
        
        % Retrieve schnitzcells and p
        s_all = myPhosphoData.(current_strain).(current_repetition).s_all
        p = myPhosphoData.(current_strain).(current_repetition).p
        
        % Repair N+1 bug
        if isfield(s_all, 'frames')
            disp('Warning: converting frames to frame_nrs.. (N+1 bug still present in your data!)')
            % Loop over all schnitzes and add the corrected frame_nrs.
            for i = 1:numel(s_all) 
                s_all(i).frame_nrs = s_all(i).frames - 1;
            end
        end

        % Make sure DJK_analyzeMu dumps data to workspace (bit hacky)
        p.dumpPlot=1;
        % Make plot
        fitTime = DJK_analyzeMu(p, s_all, 'xlim', [0 1000], 'onScreen', 0,'fitTime',[0 1150]);

        % Save dumped time & length data into myPhosphoData
        myPhosphoData.(current_strain).(current_repetition).data_time = data_time;
        myPhosphoData.(current_strain).(current_repetition).data_muField_sum = data_muField_sum;

    end
end



%% Plotting

figure(1); clf;
set(gca,'FontSize',20);
LineHandlesForLegend = [];

% Use same loop as above for plotting
% Loop over strains contained by myPhosphoData
strains = fieldnames(myPhosphoData)
for str_idx = 1:numel(strains)

    current_strain = char(strains(str_idx));
    repetitions = fieldnames(myPhosphoData.(current_strain))
    
    % Loop over repetitions done for that strain
    for rep_idx = 1:numel(repetitions)
    
        current_repetition = char(repetitions(rep_idx));        
        
        % Retrieve data 
        % (I know this is a bit redundant, but makes it more modular).
        data_time = myPhosphoData.(current_strain).(current_repetition).data_time;
        data_muField_sum = myPhosphoData.(current_strain).(current_repetition).data_muField_sum;
        
        % Plot that data
        myline = semilogy(data_time, data_muField_sum); hold on;
        set(myline,'LineStyle','-','Color',preferredcolors(str_idx,:),'LineWidth',3);

    end
    
    LineHandlesForLegend(end+1) = myline;
    
end

theLegendNames = unique(myPhosphoAuxiliary.myLegendNames); % b/c repetitions not in list of lines, get unique names.
figure(1); legend(LineHandlesForLegend, theLegendNames,'Location','best');



