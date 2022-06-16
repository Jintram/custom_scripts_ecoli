
%% User settings
FRAMENR=400;

%% Loading structure data from files

% This is my directory structure
% You can ignore this and set "fileToLoad" parameter below immediately

myRootDir = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\'; % root folder of data
positionDate = '2016-12-08_asc990_lac\'; % date data was taken, subdir in root
positionName = 'pos2crop'; % multiple colonies are measured, this identifies position, is subdir of datedir

% Load the Schnitzcells datafile
SchnitzcellsFileName = [myRootDir, positionDate, positionName,'\data\',positionName,'-Schnitz.mat']
load(SchnitzcellsFileName);
    % The schnitzcells structure is as follows:
    % schnitzcells(1) describes the life of the first cell encountered.
    % - schnitzcells(1).P gives that cell's parent (which is none, or 0, for
    % the first cell).
    % - This cell will divide, so it's two children (daughters) are given
    % by schnitzcells(1).D and schnitzcells(1).E. 
    % If e.g. schnitzcells(1).D == 3, that means schnitzcell(3) contains 
    % information about cell #1's daughter. 
    % - The frames in which this cell occured is given by
    % schnitzcells(1).frame_nrs.
    % - In each frame of the movie, cells were segmented, and where given
    % numbers without knowledge of the lineage structure. The number cells
    % were given in specific frames are stored in schnitzcells(1).cellno.
    % This is an array which corresponds to the frame numbers given by
    % schnitzcells(1).frame_nrs.
    % - Much more information about the cells can be found in the 
    % schnitzcells structure. In the code below relevant examples are
    % given.

myRootDir = [myRootDir, positionDate, positionName '\']

% To visualize with checking software:
% MW_manualcheckseg(p,'manualRange',[39:274],'overwrite',1,'assistedCorrection',1);

%% Plotting the ellips fittings
% Select which frame of the movie you'd like to see:

% Make figure
figure(1); clf; hold on;
axis equal;

% Not required - builds parameter with all coordinates of this frame
xs=[]; ys=[];

% Loop over the different schnitzcells. Each schnitzcells is an
% "individual" cell from its own birth (division) to the birth of its
% children (division). Schnitzcells are numbered by the rows of their
% struct.
for idx_s = 1:length(schnitzcells) % loop over all individuals
    
    % Look whether this individual was recorded in this frame
    frame_idx = 0;
    frame_idx = find(schnitzcells(idx_s).frame_nrs==FRAMENR);

        % If it was, then plot it to the figure
        if frame_idx ~= 0

            % TODO: there's a bug in these coordinates?!
            % x0 = schnitzcells(idx_s).rp_cenX_full(frame_idx)
            % y0 = schnitzcells(idx_s).rp_cenY_full(frame_idx)

            % Retrieve coordinates for this individual
            x0 = schnitzcells(idx_s).rp_cenX_crop(frame_idx);
            y0 = -schnitzcells(idx_s).rp_cenY_crop(frame_idx);
            
            % Get length and with of ellipse for this indivual
            ra = schnitzcells(idx_s).rp_length(frame_idx)/2;
            rb = schnitzcells(idx_s).rp_width(frame_idx)/2;

            % Get angle for individual
            ang = schnitzcells(idx_s).rp_angle(frame_idx)/360*2*pi;

            % Plot individual as ellipse using library from internet
            ellipse(ra,rb,ang,x0,y0);
            % Plot individual's center
            plot(x0,y0,'o');

            % not required - builds parameter with all coordinates of this
            % frame
            xs= [xs, x0];
            ys= [ys, y0];
            
        end

end

%% Similar to ellipses, there are also "skeletons" stored for the schnitzcells.
% These are stored separately because it's a load of data.
load([myRootDir 'data\pos2crop-skeletonData.mat']);

%% Print those skeletons for a certain frame nr.

% Set up figure
figure; clf; hold on;
axis equal;

% Go over the cells in this frame and plot the skeletons for each
for cellno = 1:numel(allExtendedSkeletons{FRAMENR})
     plot( allExtendedSkeletons{FRAMENR}{cellno}(:,2)'+allMinX{FRAMENR}(cellno),...
           -(allExtendedSkeletons{FRAMENR}{cellno}(:,1)'+allMinY{FRAMENR}(cellno)),'.');
        % allMinY and allMinX give offsets,
        % allExtendedSkeletons{FRAMENR}{cellno} hold the x,y coordinates
        % relative to this offset.
end

%% Loading microscope and segmentation images

% Note that positionName and myRootDir have to be set for this to work.
imageDir = [myRootDir,'segmentation\']; % Directory w. segm. images
fileNameSeg = [positionName, 'seg', num2str(FRAMENR), '.mat']; % Matlab name w. segmentation
load([imageDir,fileNameSeg]);

% Just to have nice colors, Lc contains the segmentation
Lc_to_show = label2rgb(Lc, @jet, 'k');%, 'shuffle');
axis tight

figure; clf;
subplottight(1,2,1);% same as subplot(1,2,1);
imshow(phsub,[]);
subplottight(1,2,2); %same as subplot(1,2,2);
imshow(Lc_to_show,[]);





