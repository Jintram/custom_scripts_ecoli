
%% Indentations 

EXPORTFOLDERFRAMES='D:\Local_Data\Dropbox\Dropbox\Filamentation_recovery_Personal\MW\figures_new\Matlab_export3\indentations\';
MICRONSPERPIXEL = 0.0431;

for imgNr = [137,147,157,167]

    %%
    h=figure(imgNr); clf; hold on;
    myImg=imread(['G:\EXPERIMENTAL_DATA_2016\2016-04-07_asc777_temperatureRecovery\pos2crop\images\pos2crop-p-2-' num2str(imgNr) '.tif']);

    myCrop=myImg(594:(594+776),284:(284+467));
    %myCrop=myImg(639:1300, 363:667); % same as fluor img
    

    imshow(myCrop,[]); hold on;
    plot([size(myCrop,2) size(myCrop,2)]-11-[0,1]./MICRONSPERPIXEL,[size(myCrop,1) size(myCrop,1)]-11,'LineWidth',4,'Color',[1 1 1]);

    saveas(h, [EXPORTFOLDERFRAMES '2016-04-07_pos2crop_fr_' num2str(imgNr) '.svg']);
    
end
