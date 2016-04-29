% renames fluo yfp images (change loop index) and copies timestamp from cfp
% image to yfp image

%*****************ADJUST************************
mydir='D:\ExperimentalDataTodo\2012-07-23\pos1\';
d_cfp=dir([mydir, '*-c-*']);
d_yfp=dir([mydir, '*-y-*']);
mysavedir='D:\ExperimentalDataTodo\2012-07-23\pos1\'; % make sure it exists!
    
%***********************************************
for runfluor=1:length(d_yfp)
    orig_yfp=imread([mydir, d_yfp(runfluor).name]);
    originfo_yfp= imfinfo([mydir, d_yfp(runfluor).name]);
    originfo_cfp= imfinfo([mydir, d_cfp(runfluor).name]);
    
    disp(['loaded ', d_yfp(runfluor).name])
    
    %change Timestamp to cfp value
    newinfo_yfp=originfo_yfp;
    newinfo_yfp.DateTime=originfo_cfp.DateTime;
    
    %change loopindex to cfp loopindex
    cfp_loopidx=d_cfp(runfluor).name(8:10);
    new_yfpname=d_yfp(runfluor).name;
    new_yfpname(8:10)=cfp_loopidx;



% save renamed yfp fluor image
% *****************************************************************************

savename=[mysavedir new_yfpname];
imwrite(orig_yfp,savename ,'TIFF','Compression','none','Description',newinfo_yfp.ImageDescription)
  disp(['saved ', new_yfpname])

end


%%

yfp552=imread([mydir, 'pos1-y-552.tif']);
yfp552info= imfinfo([mydir, 'pos1-y-552.tif']);
cfp552info= imfinfo([mydir, 'pos1-c-552.tif']);