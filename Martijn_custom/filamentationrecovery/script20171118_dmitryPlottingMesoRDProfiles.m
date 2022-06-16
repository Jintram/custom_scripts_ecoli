p1 = 'H:\EXPERIMENTAL_DATA_2017_OTHERS\dmitry_fila_recov\new results - Copy\';
p2 = 'H:\EXPERIMENTAL_DATA_2017_OTHERS\dmitry_fila_recov\old results - Copy\';

load('H:\EXPERIMENTAL_DATA_2017_OTHERS\dmitry_fila_recov\Rutger_data.mat')
% overview_all

f1 = {[p1, 'length_7.5\'],...
      [p1, 'length_10.5\'],...
      [p1, 'length_13.5\'],...
      [p1, 'length_16.5\'],...
      [p1, 'length_18\'],...
      [p1, 'length_21 350s\'],...
      [p1, 'length_24\'],...
      [p1, 'length_26 100s\'],...
      [p1, 'length_30\'],...%'length_30 370s\'      
      [p2, 'length 5_diffusion 8p2_300s\'],...
      [p2, 'length 10_diffusion 8p2_300s\'],...
      [p2, 'length 12p5_diffusion 8p2_300s\'],...
      [p2, 'length 15_diffusion 8p2_300s\'],...
      [p2, 'length 17p5_diffusion 8p2_300s\'],...
      [p2, 'length 20_diffusion 8p2_300s\'],...
      [p2, 'length 22p5_diffusion 8p2_300s\'],...
      [p2, 'length 25_diffusion 8p2_300s\'],...
      [p2, 'length 27p5_diffusion 8p2_300s\'],...
      [p2, 'length 30_diffusion 8p2_300s\'],...
      [p2, 'length 32p5_diffusion 8p2_300s\'],...
      [p2, 'length 35_diffusion 8p2_300s\'],...
      [p2, 'length 38_diffusion 8p2_300s\'],...
      }; 
      


% collect all profile values and map them to a color map of N colors.
N = 100;
A = [];
for i = 1 : length(f1)
    load([f1{i}, 'save\raw_data.mat'])
    A = [A, raw_data.density_1Dmat_avgTime_avg3_sum1];    
end
val_set = linspace(min(A), max(A), N);

%col_set = jet(N);
%col_set = makeColorMap([1 1 1],[65 148 68]./255,N);
col_set = parula(N);


h1 = figure;
hold on;

%colormap jet
%greenColorMap = makeColorMap([1 1 1],[65 148 68]./255);%,[105 189 69]./255)
%colormap(greenColorMap);
colormap(parula);

%%
allL=[];
allD={};
for i = 1 : length(f1)
    
    load([f1{i}, 'save\raw_data.mat'])

    z = raw_data.density_1Dmat_avgTime_avg3_sum1;
    % Can rescale to 1:
    % z = z/max(z);
    x = ones(size(z)) * raw_data.actual_length_um;
    y = raw_data.length_vec_um /raw_data.length_vec_um(end);
    
    for l = 1:length(z)
        
        % find color
        [val, ind] = min(abs(val_set - z(l)));
        
        % including actual z values
        % plot3(x(l), y(l), z(l), 'o', 'markerfacecolor', col_set(ind,:),  'markeredgecolor', col_set(ind,:))
        
        % flat:
        plot(x(l), y(l), 'o', 'markerfacecolor', col_set(ind,:),  'markeredgecolor', col_set(ind,:))

    end
    
    allL(end+1)=raw_data.actual_length_um;
    allD{end+1}=z;
    
end

%% 
figure(h1);
xlabel('Cell length (µm)');
ylabel('Relative location along cell');
MW_makeplotlookbetter(20);

%%
plot(overview_all(:,3), overview_all(:,2), 'o', 'markerfacecolor', [1 1 1],  'markeredgecolor', 0.9*[1 1 1])
set(gca,'Color','k')





