function [myPhosphoData, myPhosphoAuxiliary] = phospho_loaddata(myPhosphoData,myPhosphoAuxiliary,posname,posdate,groupname,ID,legendname,thesetup,thecamera);
    % Loads data in struct called myPhosphoData.
    %
    % This function requires the file 
    % 
    %
    % myPhosphoData=phospho_loadscatterdata(myPhosphoData,posname,posdate,groupname,ID,legendname);
        
    disp(['Loading ' myPhosphoAuxiliary.myRootDir,posdate,'/',posname]);
    disp(['Legend name: ', legendname]);
    
    % define p "position"
    p = DJK_initschnitz(posname,posdate,'e.coli.AMOLF','rootDir', myPhosphoAuxiliary.myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],...
        'fluor1','none','fluor2','none','fluor3','none','setup',thesetup,'camera',thecamera);
    
    % load
    % ***
    [p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);
    % ***
    % save into datastruct
    myPhosphoData.(groupname).(ID).p = p;
    
    % select which Schnitzcells to take into account (all)
    s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';
    % save into datastruct
    myPhosphoData.(groupname).(ID).s_all = s_all;
    
    
    % 
    if isfield(s_all,'muP11_all')
        mymufield = 'muP11_all'
    elseif isfield(s_all,'muP11_fitNew_all')
        mymufield = 'muP11_fitNew_all'
    else
        disp('Error!')
    end
        
    [myPhosphoData.(groupname).(ID).xvalues, myPhosphoData.(groupname).(ID).yvalues] = ...
        DJK_plot_scatterColor(p, s_all, mymufield, 'time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);

    myPhosphoAuxiliary.myLegendNames = [myPhosphoAuxiliary.myLegendNames, legendname];
    
end