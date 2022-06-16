function myPhosphoData=phospho_loadscatterdata(myPhosphoData,posname,posdate,groupname,ID,legendname);
    % loads data in struct called myPhosphoData.
    %
    % myPhosphoData=phospho_loadscatterdata(myPhosphoData,posname,posdate,groupname,ID,legendname);
        
    disp(['Loading ' posdate,'/',posname]);
    
    % define p "position"
    p = DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir',myPhosphoData.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
    
    % load
    [p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);
    % save into datastruct
    myPhosphoData.(groupname).(ID).p = p;
    
    % select which Schnitzcells to take into account (all)
    s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';
    % save into datastruct
    myPhosphoData.(groupname).(ID).s_all = s_all;
    
    
    % 
    [myPhosphoData.(groupname).(ID).xvalues, myPhosphoData.(groupname).(ID).yvalues] = ...
        DJK_plot_scatterColor(p, s_all, 'muP11_all', 'time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);

    myPhosphoData.myLegendNames = [myPhosphoData.myLegendNames, legendname];
    
end