

% size of cells is predicted to have semilog distribution
% let's check in my data

% note that not all datasets seem to have the same distribution!

%% load some data
p1 = DJK_initschnitz('pos1crop','2014-05-01','e.coli.AMOLF','rootDir','D:\MICROSCOPE_EXPERIMENTS\To_Analyze\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none','badseglistname','badseglist');
[p1,schnitzcells] = DJK_compileSchnitzImproved_3colors(p1,'quickMode',1);

% get lenghts
lengths = [schnitzcells.length_fitCoef3];
volumes = [schnitzcells.length_fitCoef3].*[schnitzcells.rp_width];

%% make histogram lengths
[myhistx,myhisty]=hist(lengths,100)

% plot normally
figure(1), clf, plot([myhistx,myhisty])

% plot semilog x
figure(2), clf, semilogx([myhistx,myhisty])

%% make histogram volumes
[myhistx,myhisty]=hist(volumes,100)

% plot normally
figure(1), clf, plot([myhistx,myhisty])

% plot semilog x
figure(2), clf, semilogx([myhistx,myhisty])

