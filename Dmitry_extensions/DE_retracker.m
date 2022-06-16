function DE_retracker(p, current_frame, problems)


%% This function allows retracking (re-connecting cells)
% 
% Shows a figure with two images: 
% Left - frame "current_frame" and right - "current_frame+1"
% After all desirebale connections have been done, it overwrites a 
% tracking file (txt). The original is saved in the folder
% "posXcrop\data\original_tracks"
%
% Controls:
% 'a'     go -20 frames
% 'f'     go +20 frames
% 's'     go -1 frame
% 'd'     go +1 frame
% 'space' to change view shcnit number/cell number;
% 'o/p'   shift cell/schnit label left/right for 2 pixels. 
%
% 'w' enables connection 1-1 cell (the same cell: 1 on the left -> 1 on the right graph)
% 'e' enables connection 1-2 cell (division: 1 on the left -> 2 on the right graph)
%
% 'q'     to quit.  The updated list is saved even with improper quitting
% (which is closing the window by pressing cross).
%
% -------------------------------------------------------------------------
% INPUTS
% problems = [Nx2]; matrix obtained with DJK_analyzeTracking(); can be
% omited (then nargin = 2).
% current_frame = int; number of frame to start with.
%
% OUTPUT 
% advanced users acan output anything they want.
% For there's no real need in output, because the result is saved
% automatically.
% I could suggest outputting new_list_num variable,it gives a matrix, which
% is a copy of the one saved in txt tracking file.
%
% Typical calls:
% DE_retracker(p, 5)
% DE_retracker(p, 100, problems)
% [new_list_num] = DE_retracker(p, 100, problems)
% -------------------------------------------------------------------------
%
% Notes.
% Colormap is hard-coded in function "plot_image";


% list of initial values of switchers and flags:
tight_flag = 0;
display_list_on_screen = 0;
flag_1 = 'cell_num';
left_screen_marked_cells = [];
right_screen_marked_cells = [];
done = 0;
shift_x = 0;



if nargin == 2
    problems = [];
    issue_with_problems = 2;
    disp(['No Problems provided.']);

elseif nargin == 3
    
    if size(problems,2) ~= 2
        disp(['The size of the "problems" matrix is wrong, should be Nx2; you have '...
            num2str(size(problems,2)) 'x' num2str(size(problems, 1))])
        disp(['Something is wrong with Problems. The program enforces Problems = []']);
        problems = [];
        issue_with_problems = 1;
    elseif (size(problems, 2) == 2) & (size(problems, 1) >= 1) 
         disp(['Using chosen problems']);
         issue_with_problems = 0;
    end
    
end

% make a folder containing all original tracking files.
if ~exist([p.tracksDir 'original_tracks\']);
    mkdir([p.tracksDir 'original_tracks\']);
end



% make a data list, from which we can uniquely correspond 
% (cell_N, frame_N) to (schnit_N);
load([p.tracksDir p.movieName '_lin.mat']);
s = schnitzcells;

for i = 1:size(s, 2)
    s(i).schnit_num = i*ones(1,length(s(i).frames));
end

cell_num_list = [s.cellno];
schnit_num_list = [s.schnit_num];
frames_list = [s.frames];

cell_frame_schnit = [cell_num_list(:), frames_list(:), schnit_num_list(:)];



% find the limits of frames that were actually tracked:
[min_frame, max_frame] = find_min_max_frame(p);


% initiate figure; can be modified to proper size.
% figure('units','normalized','outerposition',[0.25/2 0.25/2 0.75 0.75]);
% full screen:
figure('units','normalized','outerposition',[0 0 1 1]);



% plot first figure.
list_old = load_old_list(p, current_frame);
[Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);



switch issue_with_problems
    case 1
        plot_tight(1, tight_flag);
        text(5, 50, {['WARNING! Something is wrong with Problems.'] ['The program enforces Problems = [ ]']},...
            'color', 'r', 'FontWeight', 'bold', 'FontSize', 12, 'BackgroundColor', 'k')
    case 0
        plot_tight(1, tight_flag);
        text(5, 30, ['Note: Using the "Problem" matrix'],...
            'color', 'g', 'FontWeight', 'bold', 'FontSize', 12, 'BackgroundColor', 'k')
    case 2
        plot_tight(1, tight_flag);
        text(5, 30, ['Note: No Problems were provided.'],...
            'color', 'g', 'FontWeight', 'bold', 'FontSize', 12, 'BackgroundColor', 'k')
end


% loop until quit:
while ~done 
    
    % disp(['after while ' num2str(current_frame)])
    
    % get a button pressed
    waitforbuttonpress;
    key1 = double(get(gcf, 'CurrentCharacter'));
    
    
    
    switch key1
        
        % connect 1 left 1 right (no division) button 'space'
        case double('w')
            
            % collect coordinates of a click:
            [x_left, y_left] = ginput(1);
            [x_right, y_right] = ginput(1);            
            
            % feed the coordinates and the images to find the number of a
            % cell that has been clicked.
            [clicked_left, clicked_right] = find_cell(x_left,  y_left,  Lc_left,...
                x_right, y_right, Lc_right);   
            
            % feed the changes to the list:
            new_addition = [clicked_left 0 0 clicked_right];
            [new_list_txt, new_list_num, new_list_xls] = correct_old_list(list_old, new_addition, display_list_on_screen);
            
            % save the corrected list:
            save_list(p, current_frame, new_list_num, new_list_xls);
            
            color_i = rand(1,3);
            
            % collect plotting info about each selected cell:
            temp_data.N = clicked_left;
            temp_data.color = color_i;
            left_screen_marked_cells = [left_screen_marked_cells; temp_data];
            temp_data = [];
            
            temp_data.N = clicked_right;
            temp_data.color = color_i;
            right_screen_marked_cells = [right_screen_marked_cells; temp_data];
            temp_data = [];

            % plot the selected cells in both graphs:
            plot_tight(1, tight_flag)
            mark_cell(clicked_left, Lc_left, color_i)
            plot_tight(2, tight_flag)
            mark_cell(clicked_right, Lc_right, color_i)
            
            
            
            % connect 1 left 2 right (division)
        case double('e')
            
            % collect coordinates of a click:
            [x_left, y_left] = ginput(1);
            [x_right_1, y_right_1] = ginput(1);
            [x_right_2, y_right_2] = ginput(1);
            
            % feed the coordinates and the images to find the number of a
            % cell that has been clicked.
            [clicked_left, clicked_right_1, clicked_right_2] = find_cell_div(x_left,  y_left,  Lc_left,...
                x_right_1, y_right_1, x_right_2, y_right_2, Lc_right);
          
            
            % feed the changes to the list:
            new_addition = [0     clicked_left      0           clicked_right_1;...
                0           0        clicked_left   clicked_right_2];
            [new_list_txt, new_list_num, new_list_xls] = correct_old_list(list_old, new_addition, display_list_on_screen);
            
            % save the corrected list:
            save_list(p, current_frame, new_list_num, new_list_xls);

            
            color_i = rand(1,3);
            
            % collect plotting info about each selected cell:
            temp_data.N = clicked_left;
            temp_data.color = color_i;
            left_screen_marked_cells = [left_screen_marked_cells; temp_data];
            temp_data = [];
            
            temp_data.N = clicked_right_1;
            temp_data.color = color_i;
            right_screen_marked_cells = [right_screen_marked_cells; temp_data];
            temp_data = [];
            
            temp_data.N = clicked_right_2;
            temp_data.color = color_i;
            right_screen_marked_cells = [right_screen_marked_cells; temp_data];
            temp_data = [];
            
            % plot the selected cells in both graphs:
            plot_tight(1, tight_flag)
            mark_cell(clicked_left, Lc_left, color_i)
            plot_tight(2, tight_flag)
            mark_cell(clicked_right_1, Lc_right, color_i)
            mark_cell(clicked_right_2, Lc_right, color_i)
            
            
          
            % shift text left/right:
        case double('p')
            shift_x = shift_x + 2;
            [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);
            
            % when change to a new view, there is a need to plot marked
            % cells again:
            show_marked_cells(Lc_left, left_screen_marked_cells, Lc_right, right_screen_marked_cells, tight_flag);
            
         case double('o')
            shift_x = shift_x - 2;
            [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);
            
            % when change to a new view, there is a need to plot marked
            % cells again:
            show_marked_cells(Lc_left, left_screen_marked_cells, Lc_right, right_screen_marked_cells, tight_flag);
            
            
            % quit
        case double('q')
            close;
            done = 1;            
            
        
        case double(' ')
            
            % change the flag
            switch flag_1
                case 'cell_num'
                    flag_1 = 'schnit_num';
                case 'schnit_num'
                    flag_1 = 'cell_num';
            end

            [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);
            
            % when change to a new view, there is a need to plot marked
            % cells again:
            show_marked_cells(Lc_left, left_screen_marked_cells, Lc_right, right_screen_marked_cells, tight_flag);
            
            
        % next frame    
        case double('s')
            % clear the list of marked cells for plotting:
            left_screen_marked_cells = [];
            right_screen_marked_cells = [];
            
            
            current_frame = current_frame - 1;
            list_old = load_old_list(p, current_frame);
            
            [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);
            
            
        % previous frame
        case double('d')
            % clear the list of marked cells for plotting:
            left_screen_marked_cells = [];
            right_screen_marked_cells = [];
            
            
            if current_frame + 2 <= max_frame
                current_frame = current_frame + 1;
                % disp(['in the s-switch' num2str(current_frame)])
                list_old = load_old_list(p, current_frame);
                [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);
                               
            else
                [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);
                text(5, 30, ['Cannot go to higher than ' num2str(max_frame) ': out of tracked frames'],...
                    'color', 'r', 'FontWeight', 'bold', 'FontSize', 12', 'BackgroundColor', 'k')
            end
            
            
            
            % previous -20 frame
        case double('a')
            % clear the list of marked cells for plotting:
            left_screen_marked_cells = [];
            right_screen_marked_cells = [];
            
            current_frame = current_frame - 20;
            list_old = load_old_list(p, current_frame);
            
            [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);
             
            
           % next +20 frame
        case double('f')
            % clear the list of marked cells for plotting:
            left_screen_marked_cells = [];
            right_screen_marked_cells = [];
            
            if current_frame + 21 <= max_frame
                current_frame = current_frame + 20;
                list_old = load_old_list(p, current_frame);
                [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);
            else
                [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x);
                text(5, 30, ['Cannot go to higher than ' num2str(max_frame) ': out of tracked frames'],...
                    'color', 'r','FontWeight', 'bold', 'FontSize', 12, 'BackgroundColor', 'k')
            end
              
            
    end 
      
   
end

end

















% ===================================================================
% % =====================  AuXILIARY FUNCTIONs: =====================
% ===================================================================


function [Lc_left, Lc_right] = make_graph(p, current_frame, cell_frame_schnit, problems, flag_1, tight_flag, shift_x)

% LEFT FIGURE:
%load segmentations:
seg_path = [p.segmentationDir p.movieName 'seg' sprintf('%0.3d',current_frame) '.mat'];
load(seg_path);

% sometimes Lc might be not there:
if exist('Lc')
    Lc_left = Lc;
    %clear Lc;
else
    disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
    disp(['Frame ' num2str(current_frame) ' doesn''t have Lc'])
    disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
end
     



% RIGHT FIGURE:
seg_path = [p.segmentationDir p.movieName 'seg' sprintf('%0.3d',current_frame + 1) '.mat'];
load(seg_path);

% sometimes Lc might be not there:
if exist('Lc')
    Lc_right = Lc;
    clear Lc;
else
    disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
    disp(['Frame ' num2str(current_frame + 1) ' doesn''t have Lc'])
    disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
end




plot_tight(1, tight_flag)
plot_image(cell_frame_schnit, Lc_left, current_frame, problems, flag_1, shift_x)

plot_tight(2, tight_flag)
plot_image(cell_frame_schnit, Lc_right, current_frame + 1, problems, flag_1, shift_x)
end





function plot_image(cell_frame_schnit, Lc, current_frame, problems, flag_1, shift_x)

% choose colormap:

% grey (0.5)
Lc_to_show = 0.5*double(Lc > 0);

% bright grey (0.75) if a schnit is problematic:
if ~isempty(problems)
    problematic_schnitzcells = problems(problems(:,2) == current_frame, 1);
    
    % color problematic cells brighter:
    for ll = problematic_schnitzcells'
        Lc_to_show(Lc == ll) = 0.75;
    end    
else
   % disp(['There was a problem with finding a problem :) in frame ' num2str(current_frame)])
end

% other colormaps:

% label(i) -> jet(i,:)
% Lc_to_show = label2rgb(Lc, @jet, 'k');%, 'shuffle');

% label(i) -> hsv(i,:)
%Lc_to_show = label2rgb(Lc, @hsv, 'k');


cla;
imshow(Lc_to_show);
axis equal
axis tight
hold on
title(['Frame ' num2str(current_frame)]);

% find all segments (cellno) indices:
segments = unique(Lc)';
% 0 is the background segment, remove it:
segments(segments == 0)=[];



% for each segment: show its name (schnit/cell)
for cellno = segments
    
    [coord_vert_i,coord_hor_i] = find(Lc == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i);
    %text(cell_cenx_i-5, cell_ceny_i, num2str(cellno),'color','w')
    
    switch flag_1
        case 'cell_num'
            text_cellNum(cell_cenx_i, cell_ceny_i, cellno, shift_x);
        case 'schnit_num'
            % !!!!!!!!!!!!!!!!!       Extremely important    !!!!!!!!!!!!!!!!!!!!!
            % For some reason the schnit structure has frames shifted +1 compared to Lc.
            % So that Lc_frame_N actually corresponds to schnitzcell.frame(N+1)!
            text_schnitNum(cell_cenx_i, cell_ceny_i, cellno, current_frame + 1, cell_frame_schnit, shift_x)   
    end
end



% show the message schnit or cell:
switch flag_1
    case 'cell_num'
        text(1,10,'Cell numbers','Color','w','FontWeight','bold','FontSize',12);
    case 'schnit_num'
        text(1,10,'Schnit numbers','Color','g','FontWeight','bold','FontSize',12);
end


%title([num2str(size(Lc)) ', shift is ' num2str(shift_x)])

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





function [new_list_txt, new_list_num, new_list_xls] = correct_old_list(list_in, new_addition, display_list_on_screen)

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

new_list_txt = [num2str(list_in) spacer_column num2str(list_out)];

% display list_in and list_out:
if display_list_on_screen == 1    
    disp(new_list_txt);
end

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




function   text_cellNum(cell_cenx_i, cell_ceny_i, cellno, shift_x)
    text(cell_cenx_i + shift_x, cell_ceny_i, num2str(cellno),'Color','w','FontWeight','bold','FontSize',8);%'BackgroundColor',[1 1 1],
end


 

function   text_schnitNum(cell_cenx_i, cell_ceny_i, cellno, current_frame, cell_frame_schnit, shift_x)
%find a subset1, where cellno appears in the first column:    
subset_1 = cell_frame_schnit(cell_frame_schnit(:,1) == cellno, :);

%find a subset2, where current_frame appears in the second column:    
subset_2 = subset_1(subset_1(:,2) == current_frame, :);

%assign the schnit num:
if numel(subset_2) == 3
    schnit_num = subset_2(3);
elseif numel(subset_2) > 3
    disp(['Problem in relating cell to schnit. Single cell should correspond to single schnit within one frame.'])
elseif isempty(subset_2) 
    disp(['An attempt to find a schnit corresponding to cell number ' num2str(cellno) ' in frame ' num2str(current_frame) ' failed.'])
    disp(['Perhaps this frame was not processed in the tracking routine.'])
    disp(['To see the earliest tracked frame drag here "posXcrop\data\posXcrop_lin.mat" and call this: min([schnitzcells.frames])'])
end

text(cell_cenx_i + shift_x, cell_ceny_i, num2str(schnit_num),'Color','g','FontWeight','bold','FontSize',8);%'BackgroundColor',[1 1 1],

end





function  show_marked_cells(Lc_left, left_screen_marked_cells, Lc_right, right_screen_marked_cells, tight_flag)

Lc = Lc_left;
marked_cells = left_screen_marked_cells;
plot_tight(1, tight_flag)

for i = 1:length(marked_cells)
    cellno = marked_cells(i).N;
    color_i = marked_cells(i).color;
    [coord_vert_i, coord_hor_i] = find(Lc == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i);
    plot(cell_cenx_i, cell_ceny_i,'kx','MarkerSize',15,'MarkerFace','k','LineWidth',6)%,'FaceColor','w')
    plot(cell_cenx_i, cell_ceny_i,'x','color', color_i, 'MarkerSize',10,'MarkerFace', color_i,'LineWidth',3)%,'FaceColor','w') 
end




Lc = Lc_right;
marked_cells = right_screen_marked_cells;
plot_tight(2, tight_flag)

for i = 1:size(marked_cells)
    cellno = marked_cells(i).N;
    color_i = marked_cells(i).color;
    [coord_vert_i, coord_hor_i] = find(Lc == cellno);
    cell_cenx_i = mean(coord_hor_i);
    cell_ceny_i = mean(coord_vert_i);
    plot(cell_cenx_i, cell_ceny_i,'kx','MarkerSize',15,'MarkerFace','k','LineWidth',6)%,'FaceColor','w')
    plot(cell_cenx_i, cell_ceny_i,'x','color', color_i, 'MarkerSize',10,'MarkerFace', color_i,'LineWidth',3)%,'FaceColor','w') 
end

   
end


function [min_frame, max_frame] = find_min_max_frame(p)
a = [p.tracksDir p.movieName '-djk-output-*.txt'];
b = dir(a);
% last name (normally with last frame):
c1 = b(1).name;
c2 = b(end).name;

% chop off parts
c1 = c1(findstr(c1,'-output-') + 8 : end);
c2 = c2(findstr(c2,'-to-') + 4 : end);

c1(findstr(c1,'-to-') : end) = [];
c2(findstr(c2,'.txt'):end) = [];

min_frame = str2num(c1);
max_frame = str2num(c2);


end





function list_old = load_old_list(p, current_frame)

% name of a tracking file:
     name_track_file = [p.movieName '-djk-output-'...
        sprintf('%0.3d',current_frame) '-to-' sprintf('%0.3d',current_frame+1) '.txt'];   
    
    % make a backup copy if it doesnt extist yet:
    if ~(exist([p.tracksDir 'original_tracks\' name_track_file]) == 2)
        % load the list from the txt file.
        list_orig = load([p.tracksDir name_track_file]);        
        % save this list
        dlmwrite([p.tracksDir 'original_tracks\' name_track_file], list_orig, 'delimiter',' ','newline', 'pc');
    end
    
    % load the tacking file (that is "unprotected" one)
    list_old = load([p.tracksDir name_track_file]);
    
end



function plot_tight(i, tight_flag)

switch tight_flag
    
    case 1
        
        % MW: function that produces tight subplot. Code adopted from:
        % http://www.briandalessandro.com/blog/how-to-make-a-borderless-subplot-of-images-in-matlab/

        % function h = subplottight(n,m,i)
        m=2; n=1;
        [c,r] = ind2sub([m n], i);
        %ax = 
        subplot('Position', [(c-1)/m, 1-(r)/n, 1/m, 1/n]);
        %         if(nargout > 0)
        %             h = ax;
        %         end
        
    case 0
        subplot(1, 2, i);
end

end