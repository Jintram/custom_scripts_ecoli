
% some templates

%% to loop over all members of struct
r = fieldnames(myPhosphoData.s732);

for r_idx = 1:length(r)
    current_repeat = char(r(r_idx));
    myPhosphoData.s732.(current_repeat)
end

%% plotting multiple graphs

figure(4); clf; 

fr = current_branchData(bac_nr).frame_nrs;
lengths = current_branchData(bac_nr).length_fitNew;
rates = current_branchData(bac_nr).muP11_all;

[ax,hline1,hline2] = plotyy(fr,lengths,fr,rates,'semilogy','plot');

hold on; % note this has to be AFTER first plots!