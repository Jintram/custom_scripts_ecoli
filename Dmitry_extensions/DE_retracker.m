function [new_list_num] = DE_retracker(p, current_frame)

%% This function allows retracking (re-connecting cells)
% 
% Shows a figure with two images: 
% Left - frame "current_frame" and right - "current_frame+1"
% Afterr all desirebale connections have been done, it saves a new
% tracking file (txt) with addition to its original name "-edited.txt"
%
% Controls:
% 'a'     go -20 frames
% 'f'     go +20 frames
% 's'     go -1 frame
% 'd'     go +1 frame
% 'w' enables connection 1-1 cell (the same cell: 1 on the left -> 1 on the right graph)
% 'e' enables connection 1-2 cell (division: 1 on the left -> 2 on the right graph)
% 'q'     to quit.  The updated list is saved even with improper quitting
% 'space' to change view shcnit/cell
% (which is closing the window by pressing cross).





load([p.tracksDir p.movieName '_lin.mat']);
S = schnitzcells;


figure('units','normalized','outerposition',[0.25/2 0.25/2 0.75 0.75]);

[Lc_left, Lc_right] = make_graph(p, current_frame);

global_list = [];
switcher_cell_to_schnit = 1;



for l = 1:10000;
    
     name_track_file = [p.movieName '-djk-output-'...
        sprintf('%0.3d',current_frame) '-to-' sprintf('%0.3d',current_frame+1) '.txt'];
   
%     name_edited_track_file = [p.movieName '-djk-output-'...
%         sprintf('%0.3d',current_frame) '-to-' sprintf('%0.3d',current_frame+1) '-edited.txt'];
    
    name_original_track_file = [p.movieName '-djk-output-'...
        sprintf('%0.3d',current_frame) '-to-' sprintf('%0.3d',current_frame+1) '-original.txt'];   
   
    
    
    % make a backup copy of the machine-tracked list: "xxx-original.txt",
    % if doesnt extist yet:
    if ~(exist([p.tracksDir name_original_track_file])==2)
        % load the list from the txt file, that is not "original"-labeled.
        list_orig = load([p.tracksDir name_track_file]);
        % save as "xxx-original.txt"
        dlmwrite([p.tracksDir name_original_track_file], list_orig, 'delimiter',' ','newline', 'pc');
    end
        
    disp('') 
   
    % let's load the tacking file (this line could be inserted into the
    % if-statement above actually)
    list_old = load([p.tracksDir name_track_file]);
        
    

    
    
    waitforbuttonpress;
    key1 = double(get(gcf, 'CurrentCharacter'));
    
    switch key1
        
        % connect 1 left 1 right (no division) button 'space'
        case double('w')
            
            [x_left, y_left] = ginput(1);
            [x_right, y_right] = ginput(1);            
            
            [clicked_left, clicked_right] = find_cell(x_left,  y_left,  Lc_left,...
                x_right, y_right, Lc_right);   
            
            new_addition = [clicked_left 0 0 clicked_right];
            
            [new_list_txt, new_list_num, new_list_xls] = correct_old_list(list_old, new_addition);
            save_list(p, current_frame, new_list_num, new_list_xls);
            
            color_i=rand(1,3);
            
            subplottight(1,2,1)
            mark_cell(clicked_left, Lc_left, color_i)
            subplottight(1,2,2)
            mark_cell(clicked_right, Lc_right, color_i)
            
            % connect 1 left 2 right (division)
        case double('e')
            
            [x_left, y_left] = ginput(1);
            [x_right_1, y_right_1] = ginput(1);
            [x_right_2, y_right_2] = ginput(1);
            
            
          [clicked_left, clicked_right_1, clicked_right_2] = find_cell_div(x_left,  y_left,  Lc_left,...
              x_right_1, y_right_1, x_right_2, y_right_2, Lc_right);
          
          new_addition = [0     clicked_left      0           clicked_right_1;...
                          0           0        clicked_left   clicked_right_2];
                      
          [new_list_txt, new_list_num, new_list_xls] = correct_old_list(list_old, new_addition);
          save_list(p, current_frame, new_list_num, new_list_xls);
          
          color_i = rand(1,3);
          subplottight(1,2,1)
          mark_cell(clicked_left, Lc_left, color_i)
          subplottight(1,2,2)
          mark_cell(clicked_right_1, Lc_right, color_i)
          mark_cell(clicked_right_2, Lc_right, color_i)

            
        
        % quit
        case double('q')
            close;
            break
            
        
        case double(' ')
            if switcher_cell_to_schnit
                make_graph_with_shnitzcells(p, current_frame);
            else
                [Lc_left, Lc_right] = make_graph(p, current_frame);
            end
            %switch:
            switcher_cell_to_schnit = ~switcher_cell_to_schnit;
            
        % next frame    
        case double('s')
            current_frame = current_frame - 1;
            [Lc_left, Lc_right] = make_graph(p, current_frame);
            
        % previous frame
        case double('d')
            current_frame = current_frame + 1;
            [Lc_left, Lc_right] = make_graph(p, current_frame);
            
              % previous frame
        case double('a')
            current_frame = current_frame - 20;
           [Lc_left, Lc_right] = make_graph(p, current_frame);
            
                          % previous frame
        case double('f')
            current_frame = current_frame + 20;
            [Lc_left, Lc_right] = make_graph(p, current_frame);
    end 
      
   
end

end




function make_graph_with_shnitzcells(p, current_frame)

%load segmentations:
seg_path = [p.segmentationDir p.movieName 'seg' sprintf('%0.3d',current_frame) '.mat'];
load(seg_path);
Lc_left = Lc;

seg_path = [p.segmentationDir p.movieName 'seg' sprintf('%0.3d',current_frame + 1) '.mat'];
load(seg_path);
Lc_right = Lc;

subplottight(1,2,1)
plot_image(p, Lc_left, current_frame, 'schnit_num')

subplottight(1,2,2)
plot_image(p, Lc_right, current_frame + 1, 'schnit_num')

end




function [Lc_left, Lc_right] = make_graph(p, current_frame)

%load segmentations:
seg_path = [p.segmentationDir p.movieName 'seg' sprintf('%0.3d',current_frame) '.mat'];
load(seg_path);
Lc_left = Lc;

seg_path = [p.segmentationDir p.movieName 'seg' sprintf('%0.3d',current_frame + 1) '.mat'];
load(seg_path);
Lc_right = Lc;

subplottight(1,2,1)
plot_image(p, Lc_left, current_frame, 'cell_num')

subplottight(1,2,2)
plot_image(p, Lc_right, current_frame + 1, 'cell_num')
end





function plot_image(p, Lc, current_frame, flag_1)

cla;

% choose colormap:

% any label -> grey (0.5)
% Lc_to_show = 0.5*double(Lc > 0);

% label(i) -> jet(i,:)
% Lc_to_show = label2rgb(Lc, @jet, 'k');%, 'shuffle');

% label(i) -> hsv(i,:)
Lc_to_show = label2rgb(Lc, @hsv, 'k');


imshow(Lc_to_show);
axis equal
axis tight
hold on
title(['Cells. Current frame:' num2str(current_frame)]);

% find all segments (cellno) indices:
segments = unique(Lc)';
segments(segments == 0)=[];

% find each segment:
for cellno = segments
    
    [coord_vert_i,coord_hor_i] = find(Lc == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i);
    
    if strcmp(flag_1,'cell_num')
    text(cell_cenx_i - 7,cell_ceny_i,num2str(cellno),'Color','k','FontWeight','bold','FontSize',8);%'BackgroundColor',[1 1 1],
    end  
end

if strcmp(flag_1,'schnit_num')
    % load lineage:
    lin_path = [p.tracksDir p.movieName '_lin.mat'];
    load(lin_path);
    s = schnitzcells;
    %go through each schnit
    for ll = 1 : size(s,2)
        %check whether it belongs to this frame:
        if any(s(ll).frames == current_frame)
            % report its number
             text(s(ll).cenx(s(ll).frames == current_frame) - 7,...
                 s(ll).ceny(s(ll).frames == current_frame),...
                 num2str(ll), 'Color','k','FontWeight','bold','FontSize',8);
        end
    end

end

if strcmp(flag_1,'cell_num')
        text(1,10,'Cell numbers','Color','w','FontWeight','bold','FontSize',12);
elseif strcmp(flag_1,'schnit_num')
    text(1,10,'Schnit numbers','Color','w','FontWeight','bold','FontSize',12);
end

end









function [clicked_left, clicked_right] = find_cell(x_left,  y_left,  Lc_left, x_right, y_right, Lc_right)

% find all segments (cellno) indices:
segments = unique(Lc_left)';
segments(segments == 0)=[];

min_dist = 1000;
clicked_left = [];

% find each segment:
for cellno = segments
        
    [coord_vert_i, coord_hor_i] = find(Lc_left == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i); 
    
    dx = cell_cenx_i - x_left;
    dy = cell_ceny_i - y_left;
    dist = sqrt(dx^2 + dy^2);
    
    if dist < min_dist
        min_dist = dist;
        clicked_left = cellno;
    end
    
end
% =============================

% find all segments (cellno) indices:
segments = unique(Lc_right)';
segments(segments == 0)=[];

min_dist = 1000;
clicked_right = [];



% find each segment:
for cellno = segments
        
    [coord_vert_i, coord_hor_i] = find(Lc_right == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i); 
    
    dx = cell_cenx_i - x_right;
    dy = cell_ceny_i - y_right;
    dist = sqrt(dx^2 + dy^2);
    
    if dist < min_dist
        min_dist = dist;
        clicked_right = cellno;
    end
    
end


end



function [clicked_left, clicked_right_1, clicked_right_2] = find_cell_div(x_left,  y_left,  Lc_left, x_right_1, y_right_1, x_right_2, y_right_2, Lc_right)

% find all segments (cellno) indices:
segments = unique(Lc_left)';
segments(segments == 0)=[];

min_dist = 1000;
clicked_left = [];

% find each segment:
for cellno = segments
        
    [coord_vert_i, coord_hor_i] = find(Lc_left == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i); 
    
    dx = cell_cenx_i - x_left;
    dy = cell_ceny_i - y_left;
    dist = sqrt(dx^2 + dy^2);
    
    if dist < min_dist
        min_dist = dist;
        %disp('find_cell_div; clicked left')
        clicked_left = cellno;
    end
    
end
% =============================

% find all segments (cellno) indices:
segments = unique(Lc_right)';
segments(segments == 0)=[];



%----right1
min_dist = 1000;
clicked_right_1 = [];



% find each segment:
for cellno = segments
        
    [coord_vert_i, coord_hor_i] = find(Lc_right == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i); 
    
    dx = cell_cenx_i - x_right_1;
    dy = cell_ceny_i - y_right_1;
    dist = sqrt(dx^2 + dy^2);
    
    if dist < min_dist
        min_dist = dist;
        % disp('find_cell_div; clicked r1')
        clicked_right_1 = cellno;
    end
    
end


% -----right2
min_dist = 1000;
clicked_right_2 = [];

% find each segment:
for cellno = segments
        
    [coord_vert_i, coord_hor_i] = find(Lc_right == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i); 
    
    dx = cell_cenx_i - x_right_2;
    dy = cell_ceny_i - y_right_2;
    dist = sqrt(dx^2 + dy^2);
    
    if dist < min_dist
        min_dist = dist;
        % disp('find_cell_div; clicked r2')        
        clicked_right_2 = cellno;
    end
    
end


end




function mark_cell(cellno,Lc, color_i)

    [coord_vert_i, coord_hor_i] = find(Lc == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i);
    plot(cell_cenx_i, cell_ceny_i,'kx','MarkerSize',15,'MarkerFace','k','LineWidth',6)%,'FaceColor','w')
    plot(cell_cenx_i, cell_ceny_i,'x','color', color_i, 'MarkerSize',10,'MarkerFace', color_i,'LineWidth',3)%,'FaceColor','w')

end




function [new_list_txt, new_list_num, new_list_xls] = correct_old_list(list_in, new_addition)

list_out = list_in;
% new_list_xls = cell(size(list_in,1), 1 + 2*size(list_in,2));

% connectiing matrix (just a filler to see where the cahnges have been
% made):
spacer_column = repmat('      ',size(list_in,1),1);

for m = 1:size(new_addition,1)
    % origin_cell = new_addition(m, 1);
    target_cell = new_addition(m, end);
    
    % find the target cell entry in the old list:
    [ind] = find(list_in(:,4) == target_cell);
    
    % change this line in the new list:
    list_out(ind,:) = new_addition(m,:);
    
    % indicate changes in the filler:
    spacer_column(ind,:) = ['  ->  '];
    
    new_list_xls{ind,5} = '->';
end

new_list_num = list_out;

% populate xls:
% for kk = 1:numel(list_in)
%     new_list_xls{kk} = list_in(kk);
%     new_list_xls{numel(new_list_num) + size(new_list_num,1) + kk} = new_list_num(kk);
% end


% display list_in and list_out:
new_list_txt = [num2str(list_in) spacer_column num2str(list_out)];
disp(new_list_txt);

end



function save_list(p, current_frame, new_list_num, new_list_xls)

name_txt_track_file = [p.movieName '-djk-output-'...
    sprintf('%0.3d',current_frame) '-to-' sprintf('%0.3d',current_frame+1) '-edited.txt'];

dlmwrite([p.tracksDir name_txt_track_file], new_list_num, 'delimiter',' ','newline', 'pc');


% name_xls_track_file = [p.movieName '-djk-output-'...
%     sprintf('%0.3d',current_frame) '-to-' sprintf('%0.3d',current_frame+1) '-edited.xls'];
% 
% xlswrite([p.tracksDir name_xls_track_file],new_list_xls)

end



function h = subplottight(n,m,i)
    % MW: function that produces tight subplot    
    % Source: [1]
    [c,r] = ind2sub([m n], i);
    ax = subplot('Position', [(c-1)/m, 1-(r)/n, 1/m, 1/n])
    if(nargout > 0)
      h = ax;
    end
    
% [1] http://www.briandalessandro.com/blog/how-to-make-a-borderless-subplot-of-images-in-matlab/
end