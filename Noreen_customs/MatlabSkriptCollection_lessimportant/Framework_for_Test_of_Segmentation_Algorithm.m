% THIS ROUTINE PROVIDES A FRAMEWORK TO TEST DIFFERENT SEGMENTATION
% ALGORITHMS.
% GOAL: FIND A SEGMENTATION WHICH WORKS PROPERLY FOR CELLS GROWN ON RICH
% MEDIUM (-> DISPLAY SUBSTRUCTURES)
%
% STEPS:
% 1) Initiate data structure 'p' and choose segmentation algorithm.create
%    segmentation folder.
% 2) Image segmentation


% ------------------------------------------------
% ---- OPTIONS -> ADJUST ----
%SegmentationRoutine='segmentation_orig';
SegmentationRoutine='segnew';

slices=[2]; % phase contrast slices. default [1 2 3].  Try also [2]

% -----
segRange=[107 221];%107;%221;   % image which will be segmented
% segmentation parameters (only adjust with good reason)
rangeFiltSize=35; % default=35
maskMargin=20; %default=20
LoG_Smoothing=2; %default=2
minCellArea=200; %default=250
GaussianFilter=5; %default=5
minDepth=5; %defaule=5
neckDepth=2; %default=2
% ------------------------------------------------


% ------------------------------------------------
% 1) Initiate 'p' and create corresponding folders
% ------------------------------------------------
p = DJK_initschnitz('pos1crop','2013-TE-ST','e.coli.AMOLF','rootDir','D:\SegmentationAlgorithm_new_RichMedium\', ...
    'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','r','fluor2','g','fluor3','none');

stringslices=sprintf('%d',slices);
p.segmentationDir= [p.dateDir SegmentationRoutine '_slice' stringslices filesep]; % directory specific for seg. algorithm
if exist(p.segmentationDir)~=7 % if it does not exist -> create it
        [status,msg,id] = mymkdir([p.segmentationDir]);
        if status == 0
           disp(['Warning: unable to mkdir ' p.segmentationDir ' : ' msg]);
            return;
        end
end


%%
% ------------------------------------------------
% 2) Image segmentation
% ------------------------------------------------

% DJK_tracker_djk: 
if strcmp(SegmentationRoutine,'segmentation_orig')==1;
    PN_segmoviephase_3colors(p,'segRange', segRange,'slices', slices,'rangeFiltSize', rangeFiltSize,'maskMargin', ...
        maskMargin,'LoG_Smoothing', LoG_Smoothing,'minCellArea', minCellArea,'GaussianFilter', GaussianFilter, ... 
        'minDepth', minDepth,'neckDepth', neckDepth);

elseif strcmp(SegmentationRoutine,'segnew')==1;
    PN_segmoviephase_3colors(p,'segRange', segRange,'slices', slices,'rangeFiltSize', rangeFiltSize,'maskMargin', ...
        maskMargin,'LoG_Smoothing', LoG_Smoothing,'minCellArea', minCellArea,'GaussianFilter', GaussianFilter, ... 
        'minDepth', minDepth,'neckDepth', neckDepth,'medium','rich');
else
    error('Cannot find chosen segmentation algorithm.')
end