

%% note that the poisson function normally is used to describe the 
% probability on n.

n=[1:70];
lambda=2;
dt=.1;
t=10;

P=((lambda*t).^n)./(factorial(n)).*(exp(-lambda.*t));
figure; clf;
plot(n,P);
xlabel(['Number of events within time \Deltat =' num2str(t)]);
ylabel('Probability');

surface = sump(P)

%%

figure; clf; hold on;

n=70;
lambda=2;
dt=.1;
t=[0:dt:100];

for lambda=1:5
  
    %{
    % note that it is normalized for n, so I divide by surface, but not
    % sure this is allowed, since 
    P=((lambda*t).^n)./(factorial(n)).*(exp(-lambda.*t));
    surface=sum(P)*dt;
    Pnorm = P./surface;
    %}
    
    % random internet people,
    % http://www.stat.yale.edu/~pollard/Courses/241.fall97/Poisson.Proc.pdf
    % say the pdf is 
    P = exp(-lambda.*t) .* lambda.^n .* t.^(n-1) ./ factorial(n-1);
    
    %plot(t,P,'LineWidth',5-lambda+1);
    
    sum(P*dt)
    tExpected = sum(P.*t)*dt;
    disp(['For lambda = ' num2str(lambda) ', expected t value is: ' num2str(tExpected)]);
    
    %plot(t,P);

    tprime = t./tExpected;
    dtprime = dt./tExpected;
    plot(tprime,P./dtprime,'LineWidth',5-lambda+1);
    
end
    
xlabel('Waiting time');
ylabel('Probability');


%% Can we now fit our data to these functions?

DATAIDX=1;
t  = yBinsCenters;
dt = yBinsCenters(2)-yBinsCenters(1);

expandedpdf = cell2mat(arrayfun(@(x) normalizedPdf{1}(x)*ones(1,10), 1:numel(normalizedPdf{1}),'UniformOutput',0));
expandedt =   cell2mat(arrayfun(@(x) t(x)*ones(1,10)+[-dt/2:dt/9:dt/2], 1:numel(t),'UniformOutput',0));

figure; plot(expandedt,expandedpdf);

% x(1) = lambda
% x(2) = n, number of events
myFun = @(x) exp(-x(1).*t) .* x(1).^round(x(2)) .* t.^(round(x(2))-1) ./ factorial(round(x(2))-1);
errorFun = @(x) ...
    sum( (normalizedPdf{DATAIDX} - myFun(x)).^2 );

myFun = @(x) exp(-x(1).*expandedt) .* x(1).^round(x(2)) .* expandedt.^(round(x(2))-1) ./ factorial(round(x(2))-1);
errorFun = @(x) ...
    sum( (expandedpdf - myFun(x)).^2 );

%x=fmincon(errorFun,[3,60],[0,-1],1)
x=fminsearch(errorFun,[3,10])

mySolution = @(t) ((x(1)*t).^round(x(2)))./(factorial(round(x(2)))).*(exp(-x(1).*t));
mySolution = @(t) exp(-x(1).*t) .* x(1).^round(x(2)) .* t.^(round(x(2))-1) ./ factorial(round(x(2))-1);

figure; clf; hold on;
plot(t,normalizedPdf{DATAIDX});
tSolution = linspace(0,max(t),100);
dtSolution=tSolution(2)-tSolution(1);
%surface = sum(mySolution(tSolution).*dtSolution)
plot(tSolution,mySolution(tSolution));

plot(expandedt,expandedpdf)

lambda = x(1)
events = x(2)


%% 
%{
figure; hold on;

for lambda=1:10
    for n=50
        
        myFun = @(t) exp(-lambda.*t) .* lambda.^n .* t.^(n-1) ./ factorial(n-1)
        plot(tSolution,myFun(tSolution));
        
    end
end
%}




