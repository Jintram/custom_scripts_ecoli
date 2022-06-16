%%

ltarget = 2;  % target length
N = 10; % inspect N*l later

lhalf = ltarget/2;

L=ltarget:.01:N*ltarget;

nrDivisions = round(L/ltarget)
divisionLengths = L./nrDivisions;

%% plot
figure(1); clf; hold on
plot(L, divisionLengths,'LineWidth',2)
MW_makeplotlookbetter(20);
xlabel('Bacterium length (L)')
ylabel('Smallest cut length')

%% plot points of switch

discontinutities = [(ltarget+lhalf):ltarget:L(end)]

low=0;
high=max(divisionLengths)+ltarget;

for idx = 1:numel(discontinutities)
    plot([discontinutities(idx), discontinutities(idx)], [low,high],':','LineWidth',2)
end

ylim([low,high])

%% How do you expect one cell to divide?
% This is tricky, since you also still expect them to keep growing.

LENGTH=35;
STD=1;
N=1000;
colonyLengths = normrnd(LENGTH,STD,1,N);

%colonyLengths = [10];


%
observedCells = [];


division=1;
while(division)

    observedCells = [observedCells colonyLengths];
    
    % determine which ones divide
    divideIndices = find(colonyLengths > 2*ltarget);    
    
    if isempty(divideIndices)
        break
    end
    
    % divide them
    for idx = divideIndices
        
        currentLength = colonyLengths(idx);
        
        % determine # division planes
        nrDivisions = round(currentLength/ltarget);
        lresult = currentLength/nrDivisions;

        % determine where to divide
        % pick randomnumber between 1 and nrDivisions-1
        divisionN = round(.5+rand()*(nrDivisions-.5));
        
        % updates lineageLength
        cutlength = divisionN*lresult;
        colonyLengths(idx) = cutlength;
        colonyLengths(end+1) = (currentLength-cutlength);
        
    end
    
end

%observedCells

%%

figure(2); clf; hold on;
[count, x] = hist(observedCells,100);
plot(x,count,'-','LineWidth',2);
MW_makeplotlookbetter(20);
xlabel('cell length');
ylabel('count');








