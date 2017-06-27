function [t,y,normy]=script20160603_z_signalanalysis(IMGPATHSTART,FRAMERANGE,SWITCHTIME)
%% Extract time data
listexptimestr  = {};
listgainstr     = {};
listexptime     = nan(1,max(FRAMERANGE));
listcube        = {};
listdatenumber  = nan(1,max(FRAMERANGE));
for frIdx = FRAMERANGE

    % load image information    
    [exptimestr, gainstr,exptime,cube,datenumber]=...
        imsettings([IMGPATHSTART sprintf('%03d', frIdx) '.tif']);
        
    if isa(datenumber,'char')
    if strcmp(datenumber,'empty')
        datenumber=nan;
    end
    end
    
    % save image information
    listexptimestr{frIdx}   = exptimestr;
    listgainstr{frIdx}      = gainstr;
    listexptime(frIdx)      = exptime;
    listcube{frIdx}         = cube;
    listdatenumber(frIdx)   = datenumber;
    
end

% [exptimestr, gainstr,exptime,cube,datenumber]=imsettings(['H:\EXPERIMENTAL_DATA_2017\2017_04_18_switchTimeTestLong_only\switch_time_test\long_tube_test_one\pos1\pos1-g-' num2str(frIdx) '.tif'])

%%
myTimes = listdatenumber;%-min(listdatenumber(FRAMERANGE));

myTimesSinceSwitch = myTimes-SWITCHTIME;
myTimesSinceSwitchHrs = myTimesSinceSwitch*24


%% Now find the signal for each frame

signal     = nan(1,max(FRAMERANGE));
for frIdx = FRAMERANGE

    currentImg = imread([IMGPATHSTART sprintf('%03d', frIdx) '.tif']);
    signal(frIdx) = mean(currentImg(:));

end

%%
%figure;plot(myTimesSinceSwitchHrs(FRAMERANGE)*60,signal(FRAMERANGE),'ok-','LineWidth',2)

t=myTimesSinceSwitchHrs(FRAMERANGE)*60;
y=signal(FRAMERANGE);

normy = (y-min(y))./max((y-min(y)));

end


