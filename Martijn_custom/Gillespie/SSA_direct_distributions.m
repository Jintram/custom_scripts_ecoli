
% Written by Johannes Keegstra
% 2015/02

clear all

t=0;
tstop = 10; %% Specify initial and final times
x0 = [0; 20]; %% Specify initial conditions

S = [1 -1 0 0; 0 0 1 -1]; %% Specify stoichiometry

NTimes = 10;
T_Array = linspace(0,tstop,NTimes);  % times at which to plot.

prop = @(x)([10; 1*x(1); 10*x(1); 1*x(2)]);

NSims = 500;  %number of simulations;
X_Array = zeros(2,NTimes,NSims);

reverseStr='';
figure(1)
col=copper(NSims);
for isim = 1:NSims
    X_Array(:,:,isim) = SSA_direct_function(S,prop,x0,T_Array);
    
    %%  <-- this gives a cell structure.  Hit %enter to run cell.
    
    subplot(3,1,1)
    
    plot(T_Array,X_Array(1,:,isim),'-*','color',col(isim,:));  % Plot mRNA vs time
    hold on
    
    subplot(3,1,2)
    plot(T_Array,X_Array(2,:,isim),'-*','color',col(isim,:));  % Plot Protein vs time
    hold on
    
    subplot(3,1,3)
    plot(X_Array(1,:,isim),X_Array(2,:,isim),'-*','color',col(isim,:));  % Plot mRNA
    hold on
    
    % drawnow 
     msg = sprintf('Processed %d/%d', isim, NSims);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end



%%
f2 = figure(2);
BINS = [linspace(0,30,40);linspace(0,200,40)];

for iT = 1:length(T_Array)
f2 = figure(2); clf;
    for iSpec = 1:2
        subplot(2,1,iSpec)
        counts = hist(squeeze(X_Array(iSpec,iT,:)),BINS(iSpec,:));
        plot(BINS(iSpec,:),counts);
    end
    drawnow
    Movie(iT) = getframe(f2);
end

movie2avi(Movie,'tmp.avi')
    
    
