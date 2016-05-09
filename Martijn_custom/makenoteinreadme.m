function makenoteinreadme(thedirectory,scriptname);

%% write 'files added by scriptname' to scriptinfo file
fid = fopen([thedirectory 'scriptinfo.txt'], 'a');
fprintf(fid, ['\n ' datestr(datetime()) ' - Plots were added by ' scriptname]);
fclose(fid);

end