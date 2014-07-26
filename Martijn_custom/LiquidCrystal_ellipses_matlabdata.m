
%% Loading our structure from scratch

% This is my directory structure
myRootDir = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze'; % root folder of data
positionDate = '2014-05-01'; % date data was taken, subdir in root
positionName = 'pos1crop'; % multiple colonies are measured, this identifies position, is subdir of datedir

% You can edit the parameter below immediately
% Path to schnitzcells file:
fileToLoad = [myRootDir, '\', positionDate, '\',positionName,'\data\',positionName,'_lin.mat']

% Load the Schnitzcells datafile
load(fileToLoad);

% To visualize with checking software:
% MW_manualcheckseg(p,'manualRange',[39:274],'override',1,'assistedCorrection',1);

%% Plotting
frame = 120

figure(1)
clf
hold on;

xs=[];
ys=[];

for idx_s = 1:length(schnitzcells)
    
    frame_idx = 0;
    frame_idx = find(schnitzcells(idx_s).frames==frame);

        if frame_idx ~= 0

            x0 = schnitzcells(idx_s).rp_cenX_full(frame_idx)
            y0 = schnitzcells(idx_s).rp_cenY_full(frame_idx)

            x0 = schnitzcells(idx_s).rp_cenX_crop(frame_idx)
            y0 = -schnitzcells(idx_s).rp_cenY_crop(frame_idx)
            
            ra = schnitzcells(idx_s).rp_length(frame_idx)/2
            rb = schnitzcells(idx_s).rp_width(frame_idx)/2

            ang = schnitzcells(idx_s).rp_angle(frame_idx)/360*2*pi

            ellipse(ra,rb,ang,x0,y0)
            plot(x0,y0,'o');

            xs= [xs, x0];
            ys= [ys, y0];
            
        end

end

xs
ys