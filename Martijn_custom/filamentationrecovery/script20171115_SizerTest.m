
% Note that generally sizer behavior is easily recognizable beacause in
% that case cells only divide at a very narrow range of parent size; Fig 1D
% shows that this is not the case.


%DivSizes=[4:4:100];
%DivSizes=[6:6:100];
regimewidth=3;
DivSizes=[2*regimewidth:2*regimewidth:100];

TraceCollection={};
deltaCollection={};
birthSizeCollection={};
areaCollection={};
dt=.1;
theTimes=0:dt:500;
for simIdx=1:100
    
    L=40;
    birthSizes=40;
    birthTIdx=1;
    sizerMistake=0;        

    N=1;
    mu=1/60;
    deltas=[];
    areas=[];
    for tIdx=2:numel(theTimes)
        
        t=theTimes(tIdx);


        N(end+1) = N(end) + dt*(-(N(end)-1)/10 + (rand(1)-.5)*2/10);

        L(end+1) = L(end) + dt*L(end).*mu*N(end);


        %if any((L(end))>DivSizes & (L(end)-.1)<DivSizes)
            % no mistake
        %if any((L(end).*sizerMistake+.1)>DivSizes &(L(end)-.1).*sizerMistake<DivSizes) 
            % relative mistake
        if any((L(end)+sizerMistake+.1)>DivSizes & (L(end)+sizerMistake-.1)<DivSizes)
            % absolute mistake
        
            nrSites = round(L(end)/(2*regimewidth));
            theSite=floor(rand().*nrSites+1);
            frac=(theSite*2-1)./(2*nrSites);

            deltas(end+1)=L(end)-birthSizes(end);
            areas(end+1)=sum(L(birthTIdx:end).*dt);
            
            % division
            %L(end) = L(end).*frac; % no mistake
            %L(end) = L(end).*frac.*(1+(rand()-.5)/10); % proportional mistake
            L(end) = L(end).*frac+((rand()-.5).*2)*.3; % absolute mistake
            
            birthSizes(end+1)=L(end);
            birthTIdx = tIdx;
            %sizerMistake=(1+(rand()-.5)/10); % for proportional
            sizerMistake=((rand()-.5).*2)*.3; % for absolute mistake
        end

    end

    TraceCollection{end+1}=L;
    deltaCollection{end+1}=[deltas NaN];
    birthSizeCollection{end+1}=birthSizes;
    areaCollection{end+1}=[areas NaN];
end

figure(1);
plot(theTimes,TraceCollection{1},'.-')


%%
figure(2); clf; hold on;
for i=1:50
    plot(theTimes,TraceCollection{i},'.-')
end

%%
figure(3); clf; hold on; title('Added size');
scatter([birthSizeCollection{:}],[deltaCollection{:}]);

%% 
figure(4); clf; hold on; title('Area under curve');
scatter([birthSizeCollection{:}],[areaCollection{:}]);
%%

%{
theSite=[];
sizerMistake=[];
for i=1:100000
theSite(end+1)=floor(rand().*nrSites+1);
sizerMistake(end+1)=((rand()-.5).*2)*.3;
end
figure;histogram(theSite,[.5:3.5])
figure;histogram(sizerMistake)
%}


