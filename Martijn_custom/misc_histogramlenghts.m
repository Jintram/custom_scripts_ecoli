

% size of cells is predicted to have semilog distribution
% let's check in my data

% note that not all datasets seem to have the same distribution!

figure(1), clf, hold on;
figure(2), clf, hold on;

some_colors;

%% load some data
linecount=3;
switch linecount
    case 1
        p1 = DJK_initschnitz('pos1crop', '2014-05-01','e.coli.AMOLF','rootDir','D:\MICROSCOPE_EXPERIMENTS\To_Analyze\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none','badseglistname','badseglist');
    case 2
        p1 = DJK_initschnitz('pos4crop', '2014_06_18','e.coli.AMOLF','rootDir','D:\MICROSCOPE_EXPERIMENTS\To_Analyze\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none','badseglistname','badseglist');
    case 3
        p1 = DJK_initschnitz('pos2crop', '2014_06_18','e.coli.AMOLF','rootDir','D:\MICROSCOPE_EXPERIMENTS\To_Analyze\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none','badseglistname','badseglist');
end
[p1,schnitzcells] = DJK_compileSchnitzImproved_3colors(p1,'quickMode',1);

%'pos1crop', '2014-05-01'
%'pos4crop', '2014_06_18'
%'pos2crop', '2014_06_18'

% get lenghts
lengths = [schnitzcells.length_fitCoef3];
volumes = [schnitzcells.length_fitCoef3].*[schnitzcells.rp_width];

%% make histogram lengths
[myhisty,myhistx]=hist(lengths,100);
dx = myhistx(2)-myhistx(1);
myhisty_norm = myhisty./sum(myhisty.*dx);

% plot normally
figure(1), plot(myhistx,myhisty_norm,'Color',preferredcolors(linecount,:),'Linewidth',3)

set(gca,'FontSize',20);
title('Distribution of cell sizes');
xlabel('Cell length (um)');
ylabel('Count');

%% make histogram lengths on xlog scale
% plot semilog x
figure(2), semilogx(myhistx,myhisty_norm,'Color',preferredcolors(linecount,:),'Linewidth',3)

set(gca,'FontSize',20);
title('Distribution of cell sizes');
xlabel('Cell length (um)');
ylabel('Count');
%plot([min(bins_t_center), max(bins_t_center)],[0,0],'k'); % zero line


%% make histogram volumes (but that's only factor)
[myhistx,myhisty]=hist(volumes,100)

% plot normally
figure(1), plot([myhistx,myhisty])

% plot semilog x
figure(2), semilogx([myhistx,myhisty])

