function PPMs = regulationproject_abundances(strPax)

% Execute e.g. by 
% regulationproject_abundances('b4025')

databases = {'Lu', 'Lewis', 'Mancuso', 'Kim', 'Wright', 'Taniguchi'};

PPMs = NaN(1,numel(databases));

% Load dataset if not detected
if ~exist('Integrated', 'var') || ~exist('Taniguchi', 'var')
    load('U:\PROJECTS\A_NewProject\abundances-paxdb\database_all_abundances.mat');
    % WHY DOES FOLLOWING NOT WORK?
    global Integrated Kim Lewis Lu Manusco Taniguchi Wright
end

% For all separate databses
outputstr = '';
for i = 1:numel(databases)   
    eval(['hitidx = findInDB(strPax, ' databases{i} ');'])    
    if isempty(hitidx) 
        PPMs(i) = NaN;
    else
        eval(['PPMs(' num2str(i) ') = cell2mat(' databases{i} '(hitidx,3));']);
    end
    outputstr = [outputstr, num2str(PPMs(i)), '\t '];
end

sprintf(outputstr)
outputstr = strrep(outputstr, '.', ',');
sprintf(outputstr)

% For integrated one
hitidx = findInDB(strPax, Integrated);
if isempty(hitidx) 
    IntegratedPPM = NaN;
else
    IntegratedPPM = cell2mat(Integrated(hitidx,3));
end

outputstr2 = num2str(IntegratedPPM);

outputstr2 = strrep(outputstr2, '.', ',');
sprintf(['Integrated= ', outputstr2])

hitidx

end




function theIndex = findInDB(strPax,database) 

strcmp(strPax,database);
indices = strfind(database(:,2),strPax);
ix=cellfun(@isempty,indices);
indices(ix)={0}; 
indices = cell2mat(indices);

theIndex = find(indices);

end
