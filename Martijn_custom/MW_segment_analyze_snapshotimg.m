

%%%%%%%%%%%%%%%%%%%%%%%
%
% Switched to Excel sheet, this function might be redundant.
%
%%%%%%%%%%%%%%%%%%%%%%%



% 
imageToSegment = imread('F:\A_Tans1_step2_incoming_backed_up\2015-04-04\553ace\50ms\ph6.tif');


% Initialize dummy p struct.
p={}, p.movieName = 'dummy';
% Load settings into p struct. (Command taken from Excel sheet!)
p = MW_setparams_segmentation(p,'segRange', [39:274],'slices', [1 2 3],'rangeFiltSize', 35,'maskMargin', 20,'LoG_Smoothing', 2,'minCellArea', 250,'GaussianFilter', 5,'minDepth', 5,'neckDepth', 2);


% Set up input params for segmentation
% Other params to set:
saveDirectory = 'C:\Users\wehrens\Desktop\testSnapshot\';
% Code stolen from PN_segmoviephase_3colors:
inputsOfSegmentation = {'rangeFiltSize',p.rangeFiltSize,'maskMargin',p.maskMargin,...
    'useFullImage',p.useFullImage, 'LoG_Smoothing',p.LoG_Smoothing,'minCellArea',p.minCellArea,...
    'GaussianFilter',p.GaussianFilter,'minDepth',p.minDepth,...
    'neckDepth',p.neckDepth,'saveSteps',p.saveSteps,'saveDir',saveDirectory,'useFullImage',1};

% Perform the segmentation
p.useFullImage = 1;
PN_segmoviephase_3colors(p,'segRange', [6:8],'slices', [1 2 3],'rangeFiltSize', 35,'maskMargin', 20,'LoG_Smoothing', 2,'minCellArea', 250,'GaussianFilter', 5,'minDepth', 5,'neckDepth', 2);

%[phsub,LNsub,rect]= PN_segphase(imageToSegment,inputsOfSegmentation{:});

% Set

%
PN_manualcheckseg(p,'manualRange',[6:9],'overwrite',1,'assistedCorrection',1);
%MW_manualcheckseg(p,'manualRange',[39:274],'overwrite',1,'assistedCorrection',1);





%{
varargin  = {'rangeFiltSize',p.rangeFiltSize,'maskMargin',p.maskMargin,...
    'useFullImage',p.useFullImage, 'LoG_Smoothing',p.LoG_Smoothing,'minCellArea',p.minCellArea,...
    'GaussianFilter',p.GaussianFilter,'minDepth',p.minDepth,...
    'neckDepth',p.neckDepth,'saveSteps',p.saveSteps,'saveDir',saveDirectory}

varargin  = {'rangeFiltSize',p.rangeFiltSize,'maskMargin',p.maskMargin,...
    'useFullImage',p.useFullImage, 'LoG_Smoothing',p.LoG_Smoothing,'minCellArea',p.minCellArea,...
    'GaussianFilter',p.GaussianFilter,'minDepth',p.minDepth,...
    'neckDepth',p.neckDepth,'saveSteps',p.saveSteps,'saveDir',saveDirectory,'noCropFlag',true}
%}












