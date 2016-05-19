
% identify fluor colors

searchPath = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-15\pos1\';

existingNrs = [];
for framenr = 1:999
    
    fileName = ['pos1-c-' sprintf('%03d',framenr) '.tif'];

    if exist([searchPath fileName],'file')
        existingNrs(end+1) = framenr;
    end
    
end

existingNrs

mat2str(existingNrs)