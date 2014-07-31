function imagecorrel=NW_prepareImages_French_correl_average_frame(p,imageslices,frame,mysaveDir)

%n=32; %%nombre de plan dans la stack (une puissance de 2)
n=size(imageslices,3);
z=zeros([n 1]);
z(:)=(1:n);
%I=imread(['average' num2str(1,'%03d') '.tif']); %peut etre utiliser '%02d'
I=imageslices(:,:,1);
[ny nx]=size(I);
a = zeros([n 1]); %a = value of the pixel vs z position
%Im_correl = zeros([nx ny n]); size(Im_correl);
Im_correl = zeros([ny nx n]);  %vertauscht nx ny
profile =  -(z-17).*exp(-(z-17).*(z-17)/25); %%ces parametres devraient marcher pour coli


%for i=1:n
%    Im_correl(:,:,i)  = double(imread(['average' num2str(i,'%03d') '.tif']));
%end
Im_correl=double(imageslices); 

bmean = fft(profile);
bmean(1) = 0;
clear('profile','z','I');

%for x=1:nx
%    for y=1:ny
for x=1:ny
    for y=1:nx        
%for x=300:600
%    for y=600:1000
        a(:) = Im_correl(x,y,:);
        b =fft(a);
        b(1) = 0;
        cor = b.*conj(bmean);
        correal = real(ifft(cor));
        correl(:,1) = (-n/2:n/2-1);
        correl(1:n/2,2) = correal(n/2+1:n);
        correl(n/2+1:n,2) = correal(1:n/2);
        Im_correl(x,y,:) = correl(:,2);
    end
end
clear('a','b','cor','correal','correl','bmean');

Im_correl = (Im_correl - min(min(min(Im_correl))))/10;

%trouve l'image de coorelation au focus
std_focus=0;
for i=1:n
    imwrite(uint16(Im_correl(:,:,i)),[mysaveDir, p.movieName, '-correl_image-', num2str(i,'%02d'),'-', num2str(frame,'%03d'), '.tif'],'tiff','compression','none'); %blubb change!!!!
    %imwrite(uint16(Im_correl(:,:,i)),['correl_image' num2str(i,'%03d') '.tif'],'tiff','compression','none');
    std = std2(Im_correl(:,:,i));
    if (std>std_focus)
        std_focus = std;
        correl_index = i;
    end
end


%rescale l'intensite de l'image de correl au focus
alpha=double(floor(65535/(max(max(Im_correl(:,:,correl_index))))));
Im_correl(:,:,correl_index) = uint16(Im_correl(:,:,correl_index)*double(floor(65535/(max(max(Im_correl(:,:,correl_index)))))));
imwrite(uint16(Im_correl(:,:,correl_index)),[mysaveDir, p.movieName , '-correl_image_focus-', num2str(frame,'%03d'),'.tif'],'tiff','compression','none');

imagecorrel=imcomplement(Im_correl(:,:,correl_index));
imagecorrel=imagecorrel+(max2(-imagecorrel));

imwrite(uint16(imagecorrel),[mysaveDir, p.movieName , '-correl_image_focus-neg-', num2str(frame,'%03d'),'.tif'],'tiff','compression','none');

end

%%%---------------------------------------------------------------------------------------------
%%%---------------------------------------------------------------------------------------------
%ORIGINAL

%n=32; %%nombre de plan dans la stack (une puissance de 2)
%%z=zeros([n 1]);
%z(:)=(1:n);
%I=imread(['average' num2str(1,'%03d') '.tif']); %peut etre utiliser '%02d'
%[ny nx]=size(I);
%a = zeros([n 1]); %a = value of the pixel vs z position
%Im_correl = zeros([nx ny n]);
%profile =  -(z-17).*exp(-(z-17).*(z-17)/25); %%ces parametres devraient marcher pour coli%


%for i=1:n
%    Im_correl(:,:,i)  = double(imread(['average' num2str(i,'%03d') '.tif']));
%end


%bmean = fft(profile);
%bmean(1) = 0;
%clear('profile','z','I');

%for x=1:nx
%    for y=1:ny
%        a(:) = Im_correl(x,y,:);
%        b =fft(a);
% %        b(1) = 0;
%        cor = b.*conj(bmean);
%        correal = real(ifft(cor));
%        correl(:,1) = (-n/2:n/2-1);
%        correl(1:n/2,2) = correal(n/2+1:n);
%        correl(n/2+1:n,2) = correal(1:n/2);
%        Im_correl(x,y,:) = correl(:,2);
%    end
%end
%clear('a','b','cor','correal','correl','bmean');

%Im_correl = (Im_correl - min(min(min(Im_correl))))/10;

%%trouve l'image de coorelation au focus
%std_focus=0;
%for i=1:n
%    imwrite(uint16(Im_correl(:,:,i)),['correl_image' num2str(i,'%03d') '.tif'],'tiff','compression','none');
%    std = std2(Im_correl(:,:,i));
%    if (std>std_focus)
%        std_focus = std;
%        correl_index = i;
%    end
%end


%%rescale l'intensite de l'image de correl au focus
%alpha=double(floor(65535/(max(max(Im_correl(:,:,correl_index))))));
%Im_correl(:,:,correl_index) = uint16(Im_correl(:,:,correl_index)*double(floor(65535/(max(max(Im_correl(:,:,correl_index)))))));
%imwrite(uint16(Im_correl(:,:,correl_index)),'correl_image_focus.tif','tiff','compression','none');



