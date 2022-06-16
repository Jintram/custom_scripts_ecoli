% NW_addToSchntizes_lingrowth loads schnitzcells, calculates lambda of length for
% particular parts, and saves it again.
%
% *************************************************************************
% This fct is a closely adapted from DJK_addToSchntizes_mu. But here, a
% LINEAR fit fct is applied to cell lengths instead of an exponential fit.
% Linear growth rate is called: LAMBDA
% new field names contain 'lambda' instead of 'mu'. For some fields, that
% didn't have 'mu' in their name (e.g. av_fittedLength), the new field has
% and additional '_lambda' (e.g.  av_fittedLength_lambda) in its name.
%
% *************************************************************************
%
%
% DJK_compileSchnitzImproved must have been run before. For other length
% measures than rp_length, also DJK_addToSchnitzes_length has to be run
% first.
%
% First calcs lambda for complete cell cycle. Next calcs lambda for each timepoint
% with a particalur frame window (frameSizes). Calcs in 2 ways: one with
% and one without information of parent or offspring. In case without,
% lambda in beginning of cell cycle will be the same for a few timepoints. 
%
% For offspring info the following correction is applied: frameSize(end)
% end part of length schnitz is ftted and divLength calculated. Offspring
% length is added together, frameSize(end) begin part fitted and
% sum_birthLength calculated. The difference between these two is added to
% sum of offspring length, so a constant factor for each timepoints. For
% parent info a similar process is applied.
%
% Option to plot is available. Will plot data for last lengthField and last
% frameSize. Will plot mu for every 100th schnitz, for all or specific
% schnitzes only, use schnitzNum array. 
%
% OUTPUT
% 'p'   
% 'schnitzcells'      schnitzcells with length fields
%
% REQUIRED ARGUMENTS:
% 'p'
%
% OPTIONAL ARGUMENTS:
% 'schnitzName'
% 'onScreen' = 0      automatically save and close images
%              1      will ask to save and close images
%              2      will not show or save images (default)
% 'schnitzNum'        Array with schnitz numbers whose data will be plotted
%                     default: every 100th schnitz
% 'DJK_saveDir'       directory where images will be saved. 
%                     default: [p.analysisDir 'schnitzLambda\']
% 'frameSizes'        frame used for mu calc. default: [7]
%                     must be uneven!
% 'lengthFields'      length measurements used to calc mu 
%                     default: {'rp_length' 'length_fitCoef3b' 'length_fitNew'} 
%                     could still add length_fitCoef3

function [p,schnitzcells] = NW_addToSchnitzes_lingrowth(p,varargin) 

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1; functionName = 'NW_addToSchnitzes_lingrowth';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
  errorMessage = sprintf('%s\n%s',['Error width input arguments of ' functionName],['Try "help ' functionName '".']);
  error(errorMessage);
end

numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf('%s\n%s',['This input argument should be a String: ' num2str(varargin{i})],['Try "help ' functionName '".']);
      error(errorMessage);
    end
    fieldName = DJK_schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Override any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% Set default parameter values if they don't exist yet
if ~existfield(p,'schnitzName')
  p.schnitzName = [p.tracksDir,p.movieName,'-Schnitz.mat'];
end
% If onScreen, nothing is saved to disc automatically
if ~existfield(p,'onScreen')
  p.onScreen = 2;
end
if ~existfield(p,'frameSizes')
  p.frameSizes = [7]; % [3 5 7 9]
end
if ~existfield(p,'lengthFields')
  p.lengthFields = {'rp_length' 'length_fitCoef3b' 'length_fitNew'};
end
% check whether lengthFields is a cell, and not an array of char
if ~iscell(p.lengthFields)
  disp(['Warning: lengthFields is not a cell. Make sure that you use { } around text']);
end
% Save directory
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'schnitzLambda' filesep];
end
% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end
%--------------------------------------------------------------------------
  

%--------------------------------------------------------------------------
% Load Schnitzcells
%--------------------------------------------------------------------------
% Loading an existing schnitz file which contains the image-derived fields
if exist(p.schnitzName)~=2
  error(['Could not read ' p.schnitzName ' , which is required for quick mode']);
end
load(p.schnitzName);
disp(['Load from ''' p.schnitzName ''' completed...']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Override any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
if ~existfield(p,'schnitzNum')
  p.schnitzNum = [1:100:length(schnitzcells)];
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Determine which length values exist and which lambda's are gonna be calculated
%--------------------------------------------------------------------------
lengthFieldsCount = 0;
lengthFields = {};
nameFields = {};
for i = 1:length(p.lengthFields)
  f = char(p.lengthFields(i));
  if isfield(schnitzcells, f)
    lengthFieldsCount               = lengthFieldsCount + 1;
    lengthFields{lengthFieldsCount} = f;
    idx                             = strfind(f, 'length');
    nameFields{lengthFieldsCount}   = f(idx(1)+6:length(f));
    
    disp(['Will calc lambda for ' f]);
  end
end
for i = 1:length(p.frameSizes)
  disp(['Will calc lambda with frameSize ' num2str(p.frameSizes(i))]);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Loop over each schnitz
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)

  %--------------------------------------------------------------------
  % loop over each length measurement types
  %--------------------------------------------------------------------
  for num = 1:lengthFieldsCount
    lengthField   = char(lengthFields(num));
    name          = char(nameFields(num));

    %--------------------------------------------------------------------
    % get time and length of schnitz
    %--------------------------------------------------------------------
    i_time = schnitzcells(i).time/60;
    i_length = schnitzcells(i).(lengthField);
    %--------------------------------------------------------------------

    %--------------------------------------------------------------------
    % initialize with NaN (full length vectors) resp. empty vector (length
    % vector restricted to fluo time points) (init was modified 2013-11 NW)
    %--------------------------------------------------------------------
    schnitzcells(i).(['av_lambda' name])            = NaN;
    schnitzcells(i).(['av_fittedLength_lambda' name])  = NaN;
    schnitzcells(i).(['av_birthLength_lambda' name])   = NaN;
    schnitzcells(i).(['av_divLength_lambda' name])     = NaN; 
    for frameSize = p.frameSizes
      schnitzcells(i).(['lambda' num2str(frameSize) name '_all']) = [1:length(i_time)]*NaN;
      schnitzcells(i).(['lambda' num2str(frameSize) name]) =[];% [1:length(i_time)]*NaN;
      schnitzcells(i).(['lambdaP' num2str(frameSize) name '_all']) = [1:length(i_time)]*NaN;
      schnitzcells(i).(['lambdaP' num2str(frameSize) name]) =[];% [1:length(i_time)]*NaN;
    end
    % schnitzcells(i).(['av_mu' name])            = NaN;
    %schnitzcells(i).(['av_fittedLength' name])  = NaN;
    %schnitzcells(i).(['av_birthLength' name])   = NaN;
    %schnitzcells(i).(['av_divLength' name])     = NaN;
    %for frameSize = p.frameSizes
    %  schnitzcells(i).(['mu' num2str(frameSize) name '_all']) = [1:length(i_time)]*NaN;
    %  schnitzcells(i).(['mu' num2str(frameSize) name]) = [1:length(i_time)]*NaN;
    %  schnitzcells(i).(['muP' num2str(frameSize) name '_all']) = [1:length(i_time)]*NaN;
    %  schnitzcells(i).(['muP' num2str(frameSize) name]) = [1:length(i_time)]*NaN;
    %end
    %--------------------------------------------------------------------
    
    %--------------------------------------------------------------------
    % Add lambda for complete schnitz's lifetime
    %--------------------------------------------------------------------
    % Add whole cell cycle lambda when at least two measurements of length
    if length( i_time ) >= 2
      % calc the average lambda during the schnitz's lifetime 
      % length_t0 is not length at birth time but some hypothetical length
      % at t=0!
     % [schnitzcells(i).(['av_mu' name]) length_t0] = DJK_ExponentialFit( i_time, i_length); %EXP
      [schnitzcells(i).(['av_lambda' name]) length_t0] = NW_polyfit( i_time, i_length,1);   % i_length = av_lambda*i_time + length_t0 %LIN

      % calculate fitted length fitted data
      % schnitzcells(i).(['av_fittedLength' name]) = length_t0*(2.^(schnitzcells(i).(['av_mu' name])*(i_time))); EXP
      schnitzcells(i).(['av_fittedLength_lambda' name]) = schnitzcells(i).(['av_lambda' name]).*(i_time)+length_t0;

      % av_birthLength only when parent is known: depends on exact time birth
      if exist('schnitzP')
        %schnitzcells(i).(['av_birthLength' name]) = length_t0*(2.^(schnitzcells(i).(['av_mu' name])*(schnitzcells(i).birthTime/60))); EXP
        schnitzcells(i).(['av_birthLength_lambda' name]) = schnitzcells(i).(['av_lambda' name]).*(schnitzcells(i).birthTime/60) + length_t0; %LIN
      end
      
      % av_divLength only when daughter is known: depends on exact time birth
      if exist('schnitzD')
        %schnitzcells(i).(['av_divLength' name]) = length_t0*(2.^(schnitzcells(i).(['av_mu' name])*(schnitzcells(i).divTime/60))); EXP
        schnitzcells(i).(['av_divLength_lambda' name]) = schnitzcells(i).(['av_lambda' name]).*(schnitzcells(i).divTime/60) + length_t0; %LIN
      end
    end
    %--------------------------------------------------------------------
   
    %--------------------------------------------------------------------
    % loop over frameSizes
    %--------------------------------------------------------------------
    for frameSize = p.frameSizes

      %--------------------------------------------------------------------
      % get time and length of parent, offspring & sister    
      %--------------------------------------------------------------------
      clear schnitzP schnitzE schnitzD schnitzS;
      p_time = []; p_length = []; 
      de_time = []; de_length = [];
      s_time = []; s_length = [];
      age_offset = 0;

      % if schnitz has parent and sister, get its length. Correct for
      % unequal division length & measurement error at birthTime
      if (schnitzcells(i).P ~= 0) & (schnitzcells(i).S ~= 0)
        schnitzP = schnitzcells(schnitzcells(i).P); 
        schnitzS = schnitzcells(schnitzcells(i).S); 
        p_time = schnitzP.time/60;
        s_time = schnitzS.time/60;
        age_offset  = length( p_time );

        % figuring out how the division went works like this: fit lambda &
        % length for frameSize (or timeframe where both this schnitz and
        % its sister exist). Ratio is taken for fitted length at timepoint
        % in between parent's last and schnitz' first.
        data_length       = min(frameSize , min( length(s_time), length(i_time) ));
        temp_i_length     = schnitzcells(i).(lengthField)(1:data_length);
        temp_s_length     = schnitzS.(lengthField)(1:data_length);
        temp_time         = schnitzcells(i).time(1:data_length)/60;

        % if only 1 datapoint, will just use this single length datapoint as birthLength (incorrect)
        if data_length < 2 
          i_birthLength   = temp_i_length(1);
          s_birthLength   = temp_s_length(1);
        else
          % determine mu & length_t0
         % [temp_i_mu temp_i_length_t0] = DJK_ExponentialFit( temp_time, temp_i_length); EXP
         % [temp_s_mu temp_s_length_t0] = DJK_ExponentialFit( temp_time, temp_s_length); EXP
          [temp_i_lambda temp_i_length_t0] = NW_polyfit( temp_time, temp_i_length,1);   % LIN
          [temp_s_lambda temp_s_length_t0] = NW_polyfit( temp_time, temp_s_length,1);   % LIN


          % calculate fitted length data
         % i_birthLength   = temp_i_length_t0*(2^(temp_i_mu*schnitzcells(i).birthTime/60)); EXP
         % s_birthLength   = temp_s_length_t0*(2^(temp_s_mu*schnitzcells(i).birthTime/60)); EXP
         i_birthLength   = temp_i_lambda.*schnitzcells(i).birthTime/60+temp_i_length_t0 ;
         s_birthLength   = temp_s_lambda.*schnitzcells(i).birthTime/60+temp_s_length_t0 ; 
        end

        % save birthLengthRatio
        schnitzcells(i).(['birthLengthRatio_lambda' num2str(frameSize) name]) = i_birthLength/(i_birthLength+s_birthLength);
        
        % correct parent length for birthLengthRatio
        p_length          = schnitzcells(i).(['birthLengthRatio_lambda' num2str(frameSize) name])*schnitzP.(lengthField);
        
        % correct parent length for difference divLength & birthLength
        data_length         = min( length(p_length), frameSize);
        %[temp_p_mu temp_p_length_t0] = DJK_ExponentialFit(
        %p_time(end-data_length+1:end), p_length(end-data_length+1:end)); %EXP
        %p_divLength         = temp_p_length_t0*(2^(temp_p_mu*schnitzcells(i).birthTime/60)); %EXP
        [temp_p_lambda temp_p_length_t0] = NW_polyfit(p_time(end-data_length+1:end), p_length(end-data_length+1:end),1); %LIN
        p_divLength         = temp_p_lambda.*schnitzcells(i).birthTime/60+ temp_p_length_t0; %LIN
        
        data_length         = min( length(i_length), frameSize);
        %[temp_i_mu temp_i_length_t0] = DJK_ExponentialFit(
        %i_time(1:data_length), i_length(1:data_length)); EXP
        %i_birthLength       = temp_i_length_t0*(2^(temp_i_mu*schnitzcells(i).birthTime/60)); %EXP
        [temp_i_lambda temp_i_length_t0] = NW_polyfit(i_time(1:data_length), i_length(1:data_length),1); %LIN
        i_birthLength       =temp_i_lambda.*schnitzcells(i).birthTime/60 + temp_i_length_t0; %LIN

        p_length           = p_length + (i_birthLength-p_divLength);
      end

      % if schnitz has daughters, get sum of their length. Correct for
      % measurement error at birthTime
      if (schnitzcells(i).D ~= 0) & (schnitzcells(i).E ~= 0)
        schnitzD = schnitzcells(schnitzcells(i).D); 
        schnitzE = schnitzcells(schnitzcells(i).E); 
        data_length         = min( length(schnitzD.time), length(schnitzE.time) );
        de_time             = schnitzD.time(1:data_length)/60;
        de_length           = schnitzD.(lengthField)(1:data_length) + schnitzE.(lengthField)(1:data_length);

        % correct offspring length for difference divLength & sumBirthLengths
        data_length         = min( data_length, frameSize);
       % [temp_de_mu temp_de_length_t0] = DJK_ExponentialFit( de_time(1:data_length), de_length(1:data_length)); %EXP
       % de_birthLength      = temp_de_length_t0*(2^(temp_de_mu*schnitzcells(i).divTime/60)); %EXP
       [temp_de_lambda temp_de_length_t0] = NW_polyfit( de_time(1:data_length), de_length(1:data_length),1); %LIN
        de_birthLength      = temp_de_lambda.*schnitzcells(i).divTime/60+temp_de_length_t0; %LIN
       
        data_length         = min( length(i_length), frameSize);
       % [temp_i_mu temp_i_length_t0] = DJK_ExponentialFit( i_time(end-data_length+1:end), i_length(end-data_length+1:end));  % %EXP
       % i_divLength         = temp_i_length_t0*(2^(temp_i_mu*schnitzcells(i).divTime/60)); %EXP
        [temp_i_lambda temp_i_length_t0] = NW_polyfit( i_time(end-data_length+1:end), i_length(end-data_length+1:end),1); %LIN
        i_divLength         =temp_i_lambda.*schnitzcells(i).divTime/60 + temp_i_length_t0; %LIN
        de_length           = de_length + (i_divLength-de_birthLength);
      end
      %--------------------------------------------------------------------
      
      
      psde_length = cat(2,p_length , i_length, de_length);
      psde_time   = cat(2,p_time, i_time, de_time);

      %--------------------------------------------------------------------
      % Add mu with timeframe 
      %--------------------------------------------------------------------
      % Add when at least two measurements of length
      if length( i_time ) >= 2
        % loop over each timepoint of this schnitz
        for age = 1:length( i_time )
          % determine frame window without parent/offspring info
          if frameSize >= length( i_time )
            frame = [1:length(i_time)];
          else
            frame = [age - 0.5*(frameSize-1) : age + 0.5*(frameSize-1)];
            if frame(1) < 1
              frame = frame + (1 - frame(1));
            end
            if frame(end) > length( i_time )
              frame = frame - (frame(end) - length( i_time ));
            end
          end
          
          % calc lambda
         % [mu length_t0] = DJK_ExponentialFit( i_time(frame), i_length(frame)); %EXP
         % schnitzcells(i).(['mu' num2str(frameSize) name '_all'])(age) = mu; %EXP
         [lambda length_t0] = NW_polyfit( i_time(frame), i_length(frame),1); %LIN
         schnitzcells(i).(['lambda' num2str(frameSize) name '_all'])(age) = lambda; %LIN

          % determine frame window with parent/offspring info
          if frameSize >= length( psde_time )
            frame = [1:length(psde_time)];
          else
            frame = [age+age_offset - 0.5*(frameSize-1) : age+age_offset + 0.5*(frameSize-1)];
            if frame(1) < 1
              frame = frame + (1 - frame(1));
            end
            if frame(end) > length( psde_time )
              frame = frame - (frame(end) - length( psde_time ));
            end
          end
        
          % calc lambda
          %[mu length_t0] = DJK_ExponentialFit( psde_time(frame), psde_length(frame)); %EXP
          %schnitzcells(i).(['muP' num2str(frameSize) name '_all'])(age) = mu; %EXP
          [lambda length_t0] = NW_polyfit( psde_time(frame), psde_length(frame),1); %LIN
          schnitzcells(i).(['lambdaP' num2str(frameSize) name '_all'])(age) = lambda; %LIN
        end
      end
      %--------------------------------------------------------------------

      %--------------------------------------------------------------------
      % Add lambda with timeframe for fluor timepoints only 
      % remark Philippe 09-20-11 : may have to modify
      % if 2 colors are used at different timepoints 
      % Remark: only works when fluorimages are taken at same time points!
      % Takes time points of first existent fluor-color (usually fluor1).
      %--------------------------------------------------------------------
      %get existent fluor color
      if strcmp(p.fluor1,'none')==0
          fluor_frames_all=genvarname([upper(p.fluor1) '_frames_all']);
          if (i==1 && num==1), disp(['Using ' p.fluor1 ' frames to add lambda']); end
      elseif strcmp(p.fluor2,'none')==0
          fluor_frames_all=genvarname([upper(p.fluor2) '_frames_all']);
          if (i==1 && num==1), disp(['Using ' p.fluor2 ' frames to add lambda']); end
      elseif strcmp(p.fluor3,'none')==0
          fluor_frames_all=genvarname([upper(p.fluor3) '_frames_all']);
          if (i==1 && num==1), disp(['Using ' p.fluor3 ' frames to add lambda']); end
      else
          if (i==1 && num==1), disp(['You don''t have any fluor colors. Don''t know to which frames lambda should be associated. Exiting...']);end
          return
      end
      %add lambda to these frames with this fluorcolor
      FluorIndex=eval(['find(~isnan(schnitzcells(i).' fluor_frames_all '));']);
      if ~isempty(FluorIndex)
        schnitzcells(i).(['lambda' num2str(frameSize) name]) = schnitzcells(i).(['lambda' num2str(frameSize) name '_all'])(FluorIndex);
        schnitzcells(i).(['lambdaP' num2str(frameSize) name]) = schnitzcells(i).(['lambdaP' num2str(frameSize) name '_all'])(FluorIndex);
      end
      %--------------------------------------------------------------------

      %--------------------------------------------------------------------
      % If this is first run, offspring data is lacking for plot, so skipping
      %--------------------------------------------------------------------
      if p.onScreen ~= 2 & exist('schnitzD')
        if ~isfield(schnitzD, (['lambda' num2str(p.frameSizes(end)) char(nameFields(end)) '_all']))
          p.onScreen = 2;
          disp('First run, so offspring data is lacking: NOT PLOTTING!!');
        end
      end
      %--------------------------------------------------------------------

      %--------------------------------------------------------------------
      % Make plot
      %--------------------------------------------------------------------
      if p.onScreen ~= 2 & find(p.schnitzNum==i) 
        %--------------------------------------------------------------------
        % Prepare some data
        %--------------------------------------------------------------------
%         lengthField   = char(lengthFields(end));
%         name          = char(nameFields(end)); % plot last calculated data
%         frameSize     = p.frameSizes(end); % plot last frameSize

        psde_lambdaP_all  = schnitzcells(i).(['lambdaP' num2str(frameSize) name '_all']);

        if (exist('schnitzP'))
          %[mu_parent length_t0] = DJK_ExponentialFit( p_time, p_length); %EXP
          [lambda_parent length_t0] = NW_polyfit( p_time, p_length,1); %LIN
          pi_time               = cat(2, p_time, i_time(1));
          %parent_fittedLength   = length_t0*(2.^(mu_parent*pi_time)); %EXP
          parent_fittedLength   = lambda_parent*pi_time+ length_t0; %LIN
          psde_lambdaP_all          = cat(2, schnitzP.(['lambdaP' num2str(frameSize) name '_all']), psde_lambdaP_all);
        end

        if (exist('schnitzD')) & (exist('schnitzE'))
          %[mu_offspring length_t0]  = DJK_ExponentialFit( de_time, de_length); %EXP
          [lambda_offspring length_t0]  =NW_polyfit( de_time, de_length,1); %LIN
          sde_time                  = cat(2, i_time(end), de_time);
          %offSpring_fittedLength    = length_t0*(2.^(mu_offspring*sde_time)); %EXP
          offSpring_fittedLength    = lambda_offspring*sde_time + length_t0; %LIN

          data_length               = length(de_time);
          de_lambda_all                 = 0.5* ( schnitzD.(['lambda' num2str(frameSize) name '_all'])(1:data_length) + schnitzE.(['lambda' num2str(frameSize) name '_all'])(1:data_length));
          de_lambdaP_all                = 0.5* ( schnitzD.(['lambdaP' num2str(frameSize) name '_all'])(1:data_length) + schnitzE.(['lambdaP' num2str(frameSize) name '_all'])(1:data_length));
          psde_lambdaP_all              = cat(2, psde_lambdaP_all, de_lambdaP_all);
        end


        %--------------------------------------------------------------------
        % Make figure
        %--------------------------------------------------------------------    
        scrsz = get(0, 'ScreenSize');
        figureName = [p.movieDate ' ' p.movieName ' schnitz ' str4(i)];
        figureFileName = ['lambda' name '_' num2str(frameSize) '_schnitz' str4(i)];
        fig11 = figure('Position', [151 scrsz(4)-150 scrsz(3)-150 scrsz(4)-150], 'Name', figureName, 'visible','off');

        %--------------------------------------------------------------------
        % top plot: measured lengths and fitted data
        %--------------------------------------------------------------------    
        subplot(3,1,[1 2]);

        % plot fitted length
        plot(i_time, schnitzcells(i).(['av_fittedLength_lambda' name]), 'r-', 'LineWidth',2);
        % plot measured length
        hold on; plot(i_time,i_length, 'ko', 'LineWidth',2,'MarkerFaceColor',[0.8 0.8 0.8]); 

        % set & label axes
        ylim([1 10]); 
        xlim([i_time(1)-0.5 i_time(end)+0.5]);
        set(gca,'xtick',[]);
        ylabel([lengthField ' (um)'],'interpreter','none','FontWeight','bold','FontSize',10);

        %--------------------------------------------------------------------
        % top plot: parent (length divided by length ratio at division)
        %--------------------------------------------------------------------    
        if (exist('schnitzP'))
          figureName = [figureName ' (P:' str4(schnitzcells(i).P) ')'];
          % plot fitted length parent
          hold on; plot(pi_time, parent_fittedLength, 'r:', 'LineWidth',2);
          % plot measured length parent and offspring
          hold on; plot(p_time,p_length, 'ko', 'LineWidth',2);
          % display division ratio
          hold on; text(p_time(end),p_length(end)-0.25,['birthLengthRatio is ' num2str(round(100*schnitzcells(i).(['birthLengthRatio_lambda' num2str(frameSize) name]))) '%'],'interpreter','none');
        end

        %--------------------------------------------------------------------
        % top plot: sum of offspring length
        %--------------------------------------------------------------------    
        if (exist('schnitzD')) & (exist('schnitzE'))
          figureName = [figureName ' (D:' str4(schnitzcells(i).D) ' & E:' str4(schnitzcells(i).E) ')'];
          % plot fitted length offspring
          hold on; plot(sde_time, offSpring_fittedLength, 'r:', 'LineWidth',2);
          % plot measured length parent and offspring
          hold on; plot(de_time,de_length, 'ko', 'LineWidth',2);
        end

        % add title
        title(figureName,'interpreter','none','FontWeight','bold','FontSize',10);

        %--------------------------------------------------------------------
        % bottom plot: calculated lambda
        %--------------------------------------------------------------------    
        subplot(3,1,3);

        % plot calculated lambda over 5 points
        plot(i_time, schnitzcells(i).(['lambda' num2str(frameSize) name '_all']), 'k-', 'LineWidth',1);

        % plot calculated lambda over 5 points including parent
        hold on; plot(psde_time, psde_lambdaP_all, 'k:', 'LineWidth',1);

        % plot calculated lambda
        hold on; plot(i_time, ones(1,length(i_time))*schnitzcells(i).(['av_lambda' name]), 'r-', 'LineWidth',2);

        % set & label axes
        ylim([0 2.5]); 
        xlim([i_time(1)-0.5 i_time(end)+0.5]); 
        ylabel(['lambda, ' num2str(frameSize) ' frames'],'interpreter','none','FontWeight','bold','FontSize',10);
        xlabel('time (hrs)','interpreter','none','FontWeight','bold','FontSize',10); % interpreter to avoid problem with underscores

        % display lambda
        hold on; text(0.43,0.1,['av_lambda' name ' = ', DJK_setDecimalPlaces(schnitzcells(i).(['av_lambda' name]),2)],'sc','interpreter','none');

        %--------------------------------------------------------------------
        % bottom plot: parent
        %--------------------------------------------------------------------    
        if (exist('schnitzP'))
          % plot calculated lambda over 5 points
          hold on; plot(p_time, schnitzP.(['lambda' num2str(frameSize) name '_all']), 'k-', 'LineWidth',1);

          % plot calculated lambda parent
          hold on; plot(p_time, ones(1,length(p_time))*schnitzP.(['av_lambda' name]), 'r:', 'LineWidth',2);

          % display lambda
          hold on; text(0.03,0.1,['av_lambda' name ' = ', DJK_setDecimalPlaces(schnitzP.(['av_lambda' name]),2)],'sc','interpreter','none');
        end

        %--------------------------------------------------------------------
        % bottom plot: offspring
        %--------------------------------------------------------------------    
        if (exist('schnitzD')) & (exist('schnitzE'))
          % plot calculated lambda over 5 points
          hold on; plot(de_time, de_lambda_all, 'k-', 'LineWidth',1);

          % plot calculated lambda offspring
          hold on; plot(de_time, ones(1,length(de_time))*lambda_offspring, 'r:', 'LineWidth',2);

          % display lambda
          hold on; text(0.83,0.1,['av_lambda' name ' = ', DJK_setDecimalPlaces(lambda_offspring,2)],'sc','interpreter','none');
        end

        %----------------------------------------------------------------------
        % ASKING TO SAVE FIGURE
        %----------------------------------------------------------------------
        % Ask to save the figure
        if p.onScreen
          set(fig11,'visible','on');
          saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','Yes and Close','No','Yes');
          pause(0.2);
        else
          saveFigInput='Yes and Close';
        end
        if (upper(saveFigInput(1))=='Y')
            saveSameSize(fig11,'file',[p.DJK_saveDir figureFileName '.png'], 'format', 'png');
            if (strcmp(saveFigInput,'Yes and Close'))
                close(fig11);
                pause(0.2);
            end
          disp([' * Saved plot in ' figureFileName '.png']);
        end
        %----------------------------------------------------------------------
      end

    end % loop over frameSizes
    %--------------------------------------------------------------------

  end % loop over lengthField
  %--------------------------------------------------------------------
    
  %--------------------------------------------------------------------
  % Show progress
  %--------------------------------------------------------------------
  if rem(i,10) == 0, fprintf('.'); end
  if rem(i,100) == 0, fprintf(' finished schnitz %d\n', i); end
  %--------------------------------------------------------------------

end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Save new schnitzcells with pole information added
%--------------------------------------------------------------------------
save(p.schnitzName, 'schnitzcells');
disp('');
disp(['Save in ''' p.schnitzName ''' completed...']);
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function number_string = DJK_setDecimalPlaces(number, decimal_places);
number_string = sprintf(['%1.' num2str(decimal_places) 'f'], number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identical function as polyfit, but returns fit coefficient not in one
% array but as two different skalars
function [lambda,length0] = NW_polyfit(mytime, mylength,degree);
dummy=polyfit(mytime,mylength,degree);
lambda=dummy(1);
length0=dummy(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

