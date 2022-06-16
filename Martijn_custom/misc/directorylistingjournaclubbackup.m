
% Little script to give a directory listing in convenient format.

THEDIRECTORY = '\\storage01\data\AMOLF\groups\tans-group\Biophysics\2017_Journal_Club\2017_06_27_12thMeeting\';
%THEDIRECTORY = 'G:\Desktop_crap\checking\check4\';
listing = dir(THEDIRECTORY);

FinalList='';
for i=1:numel(listing)
    FinalList = [FinalList listing(i).name 10];
end

FinalList


%% Gather statistics on how many times we got articles from journals..
THEROOTDIRECTORY='\\storage01\data\AMOLF\groups\tans-group\Biophysics\2017_Journal_Club\';
listingRoot = dir(THEROOTDIRECTORY);

FinalListCells={};

% Go over multiple directories
for rr=1:numel(listingRoot)
    if ~isempty(strfind(listingRoot(rr).name,'2017'))
        currentDir = [THEROOTDIRECTORY listingRoot(rr).name '\'];
        currentListing = dir(currentDir);
        
        for i=1:numel(currentListing)
            
            % add if .pdf
            if ~isempty(strfind(currentListing(i).name,'.pdf'))
                FinalListCells{end+1} = currentListing(i).name;
            end
            
        end
        
    end
end

%{
for i=1:numel(listing)
    FinalListCells{end+1} = listing(i).name;
end
%}

% some summary stats
journalTitles={};
for i = 1:numel(FinalListCells)
    
    [match,noMatch]=regexp(FinalListCells{i},'\_','match','split');
    
    journalTitles{end+1} = noMatch{1};
end

%% brute force calculate stats
uniqueTitles = {}; uniqueTitlesCounts = [];
for i = 1:numel(journalTitles)
    
    
    currentTitle=lower(journalTitles{i});
    
    %%
    hitList=[];     
    % did we already see this journal?
    for x=1:numel(uniqueTitles)
    
        hitList(x)=strcmp(currentTitle,uniqueTitles{x});
        
    end
        
    % if already known
    if any(hitList)
        
        % increase count
        uniqueTitlesCounts(find(hitList))=uniqueTitlesCounts(find(hitList))+1;
        
    % else, add entry
    else
        
        if ~(strcmp(currentTitle,'.') | strcmp(currentTitle,'..') | strcmp(currentTitle,'suppl') | strcmp(currentTitle,'2017'))
        
            uniqueTitles{end+1} = currentTitle;
            uniqueTitlesCounts(end+1) = 1;

        end
    end
    
end

%%

figure; clf; hold on;

barh(uniqueTitlesCounts); 

set(gca, 'YTickLabel',uniqueTitles, 'YTick',1:numel(uniqueTitles))

%%

[B,I]=sort(uniqueTitlesCounts);

figure; clf; hold on;
barh(uniqueTitlesCounts(I)); 
set(gca, 'YTickLabel',uniqueTitles(I), 'YTick',1:numel(uniqueTitles(I)))





