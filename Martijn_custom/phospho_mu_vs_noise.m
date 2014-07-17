function myPhosphoData = phospho_mu_vs_noise(myPhosphoData,posname,posdate,groupname,ID,legendname)

% Config settigns ---------------------------------------------------------
MIN_NR_FRAMES = 20
% -------------------------------------------------------------------------

% loading
p = DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir',myPhosphoData.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');    
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

%s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

the_means = []; the_stds = []; the_noises = [];
for s_idx = 1:length(schnitzcells)
    if length(schnitzcells(s_idx).muP11_all) > MIN_NR_FRAMES
        % obtain values
        current_mean = mean(schnitzcells(s_idx).muP11_all);
        current_std = std(schnitzcells(s_idx).muP11_all);
        % add to array
        the_means = [the_means current_mean];
        the_stds = [the_stds current_std];
        the_noises = [the_noises current_std/current_mean];
    end
end

myPhosphoData.(groupname).(ID).the_means = the_means;
myPhosphoData.(groupname).(ID).the_stds = the_stds;
myPhosphoData.(groupname).(ID).the_noises = the_noises;

end
    

