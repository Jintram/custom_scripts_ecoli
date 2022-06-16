







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is probably redundant.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















function p = MW_setparams_segmentation(p,varargin)
% MW_setparams_segmentation(p,varargin)
% This function is to accomodate for easy taking of a snapshot segmentation.
% Use excel to call this function, which takes exactly similar arguments 
% as PN_segmoviephase_3colors.m. (In fact, this function contains exactly
% the same code as that function, but doesn't do the actual analysis.)
% 
% Example from Excel sheet command:
% MW_setparams_segmentation(p,'segRange', [39:274],'slices', [1 2 3],'rangeFiltSize', 35,'maskMargin', ,'LoG_Smoothing', 2,'minCellArea', 250,'GaussianFilter', 5,'minDepth', 5,'neckDepth', 2);

    numRequiredArgs = 1;
    if (nargin < 1) | ...
            (mod(nargin,2) == 0) | ...
            (~isSchnitzParamStruct(p))
        errorMessage = sprintf ('%s\n%s\n%s\n',...
            'Error using ==> segmoviephase:',...
            '    Invalid input arguments.',...
            '    Try "help segmoviephase".');
        error(errorMessage);
    end

    numExtraArgs = nargin - numRequiredArgs;
    if numExtraArgs > 0
        for i=1:2:(numExtraArgs-1)
            if (~isstr(varargin{i}))
                errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
                    'Error using ==> segmoviephase:',...
                    '    Invalid property ', num2str(varargin{i}), ...
                    ' is not (needs to be) a string.',...
                    '    Try "help segmoviephase".');
                error(errorMessage);
            end
            fieldName = DJK_schnitzfield(varargin{i}); % DJK 071122
            p.(fieldName) = varargin{i+1};
        end
    end

    %%%%%%%%%%%default values
    %growth medium (nb: special problems in rich medium: holes + fringed edges)
    if ~existfield(p,'medium')                                %special treatment in case of rich medium
        p.medium='normal';  
    end
    if strcmp(p.medium,'rich')==0 && strcmp(p.medium,'normal')==0
        disp('unknown kind of medium. Set to ''normal'' ');
        p.medium='normal';
    end
    %STEP A : finds a global mask and crop the image
    if ~existfield(p,'rangeFiltSize')                         %typical area for dectection of interesting features of the image
        p.rangeFiltSize = 35;
    end
    if ~existfield(p,'maskMargin')                            %additional margin of the mask : enlarge if cells missing on the edge
        p.maskMargin = 5;   % default used to be =20 but then image was imclosed not imdilated
    end
    if ~existfield(p,'useFullImage')                          % do/don't crop segmentation to ROI mask
        p.useFullImage=0;
    end
    %STEP B : find edges
    if ~existfield(p,'LoG_Smoothing')                         %smoothing amplitude of the edge detection filter
        p.LoG_Smoothing = 2;
    end
    if ~existfield(p,'minCellArea')                           %minimum cell area (objects smaller than that will be erased)
        p.minCellArea = 250;
    end
    % if rich medium, minimal cell area should be set to a small value (100).
    % enforce it
    if strcmp(p.medium,'rich')==1 & p.minCellArea>100
        p.minCellArea=100;
        disp('...')
        disp(['For better segmentation in rich medium, minCellArea should be small. It was reset to 100. ' ...
            ' A larger default value can only be written directly into the function.'])
        disp('...')
    end

    %STEP C : prepare seeds for watershedding
    if ~existfield(p,'GaussianFilter')                        %smoothing of the original image to find local minima within cells
        p.GaussianFilter = 5;
    end
    if ~existfield(p,'minDepth')                              %minimum accepted depth for a local minimum
        p.minDepth = 5;
    end
    %STEP E: treatment of long cells
    if ~existfield(p,'neckDepth')                             %minimum neck width to cut a too long cell
        p.neckDepth = 2;
    end
    %saving images
    if ~existfield(p,'saveSteps')                             %indicate if you want to save intermediate images
        p.saveSteps = true;
    end
    if ~existfield(p,'PN_saveDir') & p.saveSteps              %subfolder of p.segmentationDir where image treatment steps are saved
        p.PN_saveDir = ['param' '_Marg' num2str(p.maskMargin) '_LoG' num2str(p.LoG_Smoothing) '_Area' num2str(p.minCellArea) '_Depth' num2str(p.minDepth) '_Neck' num2str(p.neckDepth) filesep];
        [message errmsg] = sprintf(['Image saving folder automatically generated: ' p.PN_saveDir]);
        disp(message);
    end

    %acquisition technique: 'phasecontrast' (conventional) or 'brightfield' (French
    %style) possible
    if ~existfield(p,'method')                                %standard technique
        p.method='phasecontrast';  
    end
    if strcmp(p.method,'phasecontrast')==0 && strcmp(p.method,'brightfield')==0
        disp('unknown kind of acquisition method. Set to ''phasecontrast'' ');
        p.method='phasecontrast';
    end
    %saveDir for averaged BrightField Images
    if ~existfield(p,'NW_brightfieldsaveDir') %subfolder of p.imageDir
        p.NW_brightfieldsaveDir = ['averagedimages'];
    end
    %saveDir for correlated BrightField Images
    if ~existfield(p,'NW_correlatedsaveDir') %subfolder of p.imageDir
        p.NW_correlatedsaveDir = ['correlatedimages'];
    end
    %quickMode
    if ~existfield(p,'quickMode') %extract already averaged images
        p.quickMode = 0;
    end

end

