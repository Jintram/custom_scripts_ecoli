

%% General
%myim=imread('20140501.tif');         
%myim=imread('20140501_t2.tif');       
myim=imread('20140501_t3.tif');       
%myim=imread('dirty_bacteria_t1.tif'); 
%myim=imread('dirty_bacteria_t2.tif'); 

myim=normalize(double(myim)); % always normalize

myiminv=imcomplement(myim);
myimbw=im2bw(myim);

%% 

myimin = im2bw(myiminv,.7);

mySE = strel('disk', 1);
improc = imerode(myimin, mySE);

improc = bwconncomp(improc); % 
    % this now is an instance with multiple member variables. 
    % importantly, it contains lists of pixels belonging to one conn. comp.
improc = labelmatrix(improc); % labels connected area with +1 lbl
improc = label2rgb(improc, @jet, 'c', 'shuffle'); % gives it actual colors

figure(1);
subplottight(1,3,1);

% manually selecting regions
%imman = bwselect; % works on current axes

imshow(myimin,[]);
subplottight(1,3,2);
imshow(improc,[]);


subplottight(1,3,3);
imshow(imman,[]);



