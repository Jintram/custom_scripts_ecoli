% MUTUAL INFORMATION based on binning -> overestimate (at least heuristically in my case).
% see also Kraskov Grasberger PRE 2004 and even more applicable for my
% case: Kernel density estimators (KDE)

% **** ADJUST ***
% input vectors
rho=0.5;%corrcoeff
invec1=rho*unimu+sqrt(1-rho^2)*uniy;%rho*muRnd+sqrt(1-rho^2)*yrndpart;
invec2=unimu; % careful! role of '3' in explained variance script
numbinsvec=[1:10 ];%unique(ceil(length(invec1)./[1:200 210:10:2000]));
% ***************

%vector with I_12 (mutual information) dependent on the # bins
Infovec=zeros(size(numbinsvec));

for counter=1:length(numbinsvec)
    numbins=numbinsvec(counter);

% use bins equidistant
minvec1=min(invec1);
maxvec1=max(invec1);
binsize1=(maxvec1-minvec1)/numbins;
binrange1=minvec1:binsize1:maxvec1; %not sure if it refers to left edge or what else

minvec2=min(invec2);
maxvec2=max(invec2);
binsize2=(maxvec2-minvec2)/numbins;
binrange2=minvec2:binsize2:maxvec2; %not sure if it refers to left edge or what else


% perform binning
binnedvec1=zeros(size(binrange1)); % they are not used...
binnedvec2=zeros(size(binrange2));
binnedvecboth=zeros(length(binrange1),length(binrange2)); %row: vec1 cchanges, column: vec2 changes
prob1=zeros(length(binnedvec1)-1,1);
prob2=zeros(length(binnedvec2)-1,1);
probboth=zeros(length(binnedvecboth)-1,length(binnedvecboth)-1);


for vec1run=1:length(binrange1)-1
    idx1=find(invec1>=binrange1(vec1run)&invec1<binrange1(vec1run+1));
    prob1(vec1run)=length(idx1)/length(invec1);
end
   
for vec2run=1:length(binrange2)-1
    idx2=find(invec2>=binrange2(vec2run)&invec2<binrange2(vec2run+1));
    prob2(vec2run)=length(idx2)/length(invec2);
end

for vec1run=1:length(binrange1)-1
    for vec2run=1:length(binrange2)-1
        idx1=find(invec1>=binrange1(vec1run)&invec1<binrange1(vec1run+1));
        idx2=find(invec2>=binrange2(vec2run)&invec2<binrange2(vec2run+1));
        overlap=ismember(idx1,idx2);
        idxjoint=idx1(overlap);
        probboth(vec1run,vec2run)=length(idxjoint)/length(invec1);
    end
end

% mutual information
% I_12=SUM_x SUM_y [probboth*log(probboth/(prob1*prob2)]
I_12=0;
%sum over all entries but ignore the empty ones
for vecrun1=1:length(prob1)
    for vecrun2=1:length(prob2)
        if probboth(vecrun1,vecrun2)~=0
            I_12=I_12+probboth(vecrun1,vecrun2)*log(probboth(vecrun1,vecrun2)/(prob1(vecrun1)*prob2(vecrun2)));
        end
    end
   
end
Infovec(counter)=I_12;

end

%figure
clf
plot(numbinsvec,Infovec,'.-')
xlabel('#bins')
ylabel('mutual information')