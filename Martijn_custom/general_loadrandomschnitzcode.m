


p1 = DJK_initschnitz('pos1crop','2014-05-01','e.coli.AMOLF','rootDir','D:\MICROSCOPE_EXPERIMENTS\To_Analyze\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none','badseglistname','badseglist');
[p1,schnitzcells] = DJK_compileSchnitzImproved_3colors(p1,'quickMode',1);