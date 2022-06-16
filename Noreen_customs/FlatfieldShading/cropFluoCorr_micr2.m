% crops the fluo correction images to the ROI (important for microscope 2
% where several ROIs are used)
% Fluo corr image has to be saved manually!

% ** ADJUST ****************************
% fluo corr .mat file
myfile='D:\ExperimentalDataTodo\ShadingFlatfield\Correction_Micr2_2048x2048_mCherry10ms.mat';

origROI=[0,0,2048,2048];  % program will fail if the orig images don't have this size!
newROI=[328,328,1392,1392];

% !!!!! note the notation: col and row switched!!!!!
%    (startcol(!),startrow(!),width,height)

%default ROI settings micoscope2
% 1 == ROI(0,0,2048,2048);
% 2 == ROI(328,504,1392,1040);
% 3 == ROI(960,960,128,128);
% 4 == ROI(896,896,256,256);
% 5 == ROI(768,768,512,512);
% 6 == ROI(512,512, 1024,1024);
% 7 == ROI(328,328,1392,1392);

% **************************************

% get ROI coordinates
origROIstartrow=origROI(2)+1; %(!)
origROIstartcol=origROI(1)+1; %(!)
origROIrows=origROI(4); %(!)
origROIcols=origROI(3); %(!)
newROIstartrow=newROI(2)+1; %(!)
newROIstartcol=newROI(1)+1; %(!)
newROIrows=newROI(4); %(!)
newROIcols=newROI(3); %(!)

shiftstartrow=newROIstartrow-origROIstartrow; % the +1 shouldn't matter for the shifts
shiftstartcol=newROIstartcol-origROIstartcol;
if shiftstartrow<0 | shiftstartcol<0
    error('cannot shift to larger ROI!')
end
newmaxrow=shiftstartrow+newROIrows;
newmaxcol=shiftstartcol+newROIcols;

% find variables in fluocorr .mat file
fluocorrstrings=who('-file',myfile);
load(myfile);

% check whether all matrices have the correct start size
for i=1:length(fluocorrstrings)
    eval(['roicorrect=sum(size(' fluocorrstrings{i} ')==[origROIrows origROIcols]);'])
    if roicorrect~=2
        error('ROI/size of input matrices is wrong.')
    end
end
clear roicorrect

% crop fluocorr images
% ** This should be general to the case that the original ROI is not the
% full field but is not tested yet (only tested for origROI=fullfield=[0,0,2048,2048]!! **
for i=1:length(fluocorrstrings)
    eval([fluocorrstrings{i} '= ' fluocorrstrings{i} '(shiftstartrow+1:newmaxrow,shiftstartcol+1:newmaxcol);'])
    % the +1 is important for the startrow/col(!) - scenario: (0,0) as
    % start -> not possible
end



