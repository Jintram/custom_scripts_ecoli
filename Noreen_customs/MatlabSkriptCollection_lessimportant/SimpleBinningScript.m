varx=dC;%xx;%dC;   % bin according to this variable
vary=mu*100;%yy;%mu;

bins=[-0.8:0.1:1]; % very specfiic
binmeans=[-0.75:0.1:0.95];
ymeans=zeros(size(binmeans));

for i=1:length(bins)-1
    idx=find(varx>bins(i) &  varx<bins(i+1));
    ymeans(i)=mean(vary(idx));
end
    
plot(binmeans,ymeans,'.')