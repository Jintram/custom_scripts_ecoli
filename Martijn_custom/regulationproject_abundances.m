function PPMs = regulationproject_abundances(strPax)
% Return ppm value for enzyme measured by different databases.
% Input is id from PaxDB database. (Look for gene name in database%
% to obtain that ID.
% 
% Execute e.g. by 
% regulationproject_abundances('b4025')

% The available databases (excl. Integrated one)
databases = {'Lu', 'Lewis', 'Mancuso', 'Kim', 'Wright', 'Taniguchi'};

% Empty output array
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
    % Look for ID
    eval(['hitidx = findInDB(strPax, ' databases{i} ');'])    
    if isempty(hitidx) 
        % If not found, flag by NaN
        PPMs(i) = NaN; 
    else
        % otherwise add number
        eval(['PPMs(' num2str(i) ') = cell2mat(' databases{i} '(hitidx,3));']); 
    end
    outputstr = [outputstr, num2str(PPMs(i)), '\t '];
end

% Print to user
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

% index in matlab table
hitidx 

end




function theIndex = findInDB(strPax,database) 
% Return index of line that contains gene with id strPax in database

%strcmp(strPax,database);

% Look for strPax id in database
indices = strfind(database(:,2),strPax);
% Fill empty cells with zeros and convert to matrix
ix=cellfun(@isempty,indices);
indices(ix)={0}; 
indices = cell2mat(indices);

% Spit out index
theIndex = find(indices);

end
