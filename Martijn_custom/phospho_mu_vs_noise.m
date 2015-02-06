function myPhosphoData = phospho_mu_vs_noise(myPhosphoData,groupname,ID)
% function myPhosphoData = phospho_mu_vs_noise(myPhosphoData,myPhosphoAuxiliary,posname,posdate,groupname,ID,legendname)
% 
% Calculates the_means, the_stds and the_noises, and adds that to
% myPhosphoData.

% Config settigns ---------------------------------------------------------
MIN_NR_FRAMES = 20
% -------------------------------------------------------------------------

% loading
%{
p = DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir',myPhosphoAuxiliary.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');    
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);
%}
p = myPhosphoData.(groupname).(ID).p;
schnitzcells = myPhosphoData.(groupname).(ID).s_all;
theMuField = myPhosphoData.(groupname).(ID).muFieldName;

%s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

the_means = []; the_stds = []; the_noises = [];
for s_idx = 1:length(schnitzcells)
    if length(schnitzcells(s_idx).(theMuField)) > MIN_NR_FRAMES
        % obtain values
        current_mean = mean(schnitzcells(s_idx).(theMuField));
        current_std = std(schnitzcells(s_idx).(theMuField));
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
    

