

mylist = dir(p.movieDir)

for i = 1:numel(mylist)
    
    currentFileName = mylist(i).name;
    
    % if this file is a microscope image, move to proper directory
    if any(strfind(currentFileName,ourSettings.positionName)==1)
        frompath    = [p.movieDir currentFileName];
        topath      = [p.imageDir currentFileName];
        movefile(frompath, topath);
        disp(['Moved file, destination: ' topath]);
    end
end
