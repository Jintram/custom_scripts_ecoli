% NW_addToSchnitzes_smoothedFields_Rates
%
 %*****************************%
% TODO: SMOOTH LETZTEN EINTRAG BEI SMOOTHED FLUORRATES!!!!! TOCHTERRATES
% EINFACH UNGESMOOTHED UEBERNOMMEN!!!!
% ***************************************
% adds several fields:
%
% - length_fitNew_smooth3: (3=default, #smooth points) smoothed version of
% length_fitNew. smoothing is
% done with "smooth" and option "lowess" (weighed 1st order polynomial
% fitting, cf Locke2011,Science,SOM). Default #points used that contributed
% to smoothing. (typically ~10% of points. default: 3pts)
%
% - Y5_mean_smooth (same procedure as above, typically 50-70% of points
% contribute (strong smoothing! default: 3pts).
% The color ('y') is specified by colorfied (input)
% 
% - dY5_mean_len_dt: fluorescence production rate. total fluorescence is
% proxied by concentration*length (width doesn't change), i.e.
% Y5_mean*length_fitNew, rate calculated by division by time increment (cf
% DJK_addToSchnitzes_fluorRate_phase)
%
% - dY5_mean_len_smooth3_dt: (3=default,#smoothing points) as above, but uses smoothed data. 
%   **** Currently also adds last data field which is the difference of fluorescence of daughter
% fields and last field of parent. Maybe this destroys smoothing effect. Check in code which
% version is currently used.*****
%
% - Y5_mean_len, Y5_mean_len_smooth3: product of Y5_mean and
% length_fitNew resp. product of Y5_mean_smooth3 and length_fitNew_smooth3.
% Used for the above rate calculations. ***** First smoothing, then taking the
% product! *****
%
% **********
% RANDOM FLUORESCENCE DATA CAN BE ADDED
% **********
%
% ******
% This function has to be run for each color seperately. The smoothed
% length field is overwritten each time.
%
% Sanders suggestion: STORES AT Y_TIME !!!!
%
% OUTPUT
% 'schnitzcells'      schnitzcells
%
% REQUIRED ARGUMENTS:
% 'schnitzcells'      schnitzcells
% 'colorfield'        which fluor color
%
function [schnitzcells] = NW_addToSchnitzes_smoothedFields_Rates(p,colorfield) 

% ********* ADJUST ********************
smoothNrPts=3; % default: average over 3 points.
ADD_RND_DATA=0;
% *************************************

smoothNrPts=smoothNrPts+1; % otherwise "lowess" doesn't work.. bias seems to be necessary.


schnitzname = [p.tracksDir,p.movieName,'-Schnitz.mat'];
load(schnitzname);

% generate field names

Fluor_mean = [upper(colorfield) '5_mean']; %conc (input)
mylength=['length_fitNew']; % length (input)
Fluor_mean_smooth = [ upper(colorfield) '5_mean_smooth' num2str(smoothNrPts-1)]; % conc (smooth)
length_smooth = ['length_fitNew_smooth' num2str(smoothNrPts-1)]; % length (smooth)

Fluor_mean_length = [upper(colorfield) '5_len']; %conc*length
Fluor_mean_length_smooth = [upper(colorfield) '5_len_smooth' num2str(smoothNrPts-1)]; %conc*length (smooth)
dFluor_mean_length_dt = ['d' upper(colorfield) '5_len_dt']; %conc*length / Delta(T)
dFluor_mean_length_smooth_dt = ['d' upper(colorfield) '5_len_smooth' num2str(smoothNrPts-1) '_dt']; % conc*length/ Delta(T) (smooth)
fluortime = [upper(colorfield) '_time'];

if ADD_RND_DATA
    %random fields
    randLength=['rand_len']; % rand, length of #frames
    randFluor=['rand_fluor']; % rand, length of #fluo frames
    randLength_smooth=['rand_len_smooth' num2str(smoothNrPts-1)];
    randFluor_smooth=['rand_fluor_smooth' num2str(smoothNrPts-1)];
    randLengthFluor_smooth=['rand_fluor_len_smooth' num2str(smoothNrPts-1)];
    drandLengthFluor_smooth_dt=['d_rand_fluor_len_smooth' num2str(smoothNrPts-1) '_dt'];

    %half random field (real smoothed length, random fluo)
    randrealLengthFluor_smooth=['rand_fluor_real_len_smooth' num2str(smoothNrPts-1)];
    drandrealLengthFluor_smooth_dt=['d_rand_fluor_real_len_smooth' num2str(smoothNrPts-1) '_dt'];
end

% ***
%d_sum = ['d' upper(colorfield) '5_sum'];
%d_sum_dt = ['d' upper(colorfield) '5_sum_dt'];
%sum_all = [upper(colorfield) '5_sum_all'];

%phase_atfluor= ['phase_at' upper(colorfield)];
% ***
%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
ls = length(schnitzcells);
for ii = 1:ls
  schnitzcells(ii).(Fluor_mean_length)               = [];
  schnitzcells(ii).(Fluor_mean_length_smooth)        = [];
  schnitzcells(ii).(dFluor_mean_length_dt)            = [];
  schnitzcells(ii).(dFluor_mean_length_smooth_dt)     = [];
  schnitzcells(ii).(Fluor_mean_smooth)                = [];
  schnitzcells(ii).(length_smooth)                    = [];
  
  if ADD_RND_DATA

      % clear Random length and Fluo data
      schnitzcells(ii).(randLength_smooth)=[];
      schnitzcells(ii).(randFluor_smooth)=[];
      schnitzcells(ii).(randLengthFluor_smooth)=[];
      schnitzcells(ii).(drandLengthFluor_smooth_dt)=[];

      % clear half random
      schnitzcells(ii).(randrealLengthFluor_smooth)=[];
      schnitzcells(ii).(drandrealLengthFluor_smooth_dt)=[];
  end

  % ***  
%schnitzcells(ii).(d_sum)         = [];
%schnitzcells(ii).(d_sum_dt)      = [];
% **

end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES AND ADD SMOOTHED + TOTAL FLUO FIELDS + RANDOMFIELDS
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  s = schnitzcells(i);
  
  %------------------------------------------------------------------------
  % CREATE SMOOTHED DATA FIELDS
  %------------------------------------------------------------------------
  s.(Fluor_mean_smooth)=smooth(s.(Fluor_mean),smoothNrPts,'lowess')';
  s.(length_smooth)=smooth(s.(mylength),smoothNrPts,'lowess')';  
  if ADD_RND_DATA
      s.(randFluor_smooth)=smooth(s.(randFluor),smoothNrPts,'lowess')';
      s.(randLength_smooth)=smooth(s.(randLength),smoothNrPts,'lowess')';
  end
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % CALCULATE NEW PROXIES FOR TOTAL FLUOR CONCENTRATION
  %------------------------------------------------------------------------
   if length(s.(fluortime)) > 0
        for f = 1:length(s.(fluortime))
          age                 = find(s.time == s.(fluortime)(f));
          s.(Fluor_mean_length)(f)=s.(Fluor_mean)(f)*s.(mylength)(age);
          s.(Fluor_mean_length_smooth)(f)=s.(Fluor_mean_smooth)(f)*s.(length_smooth)(age);
          if ADD_RND_DATA
              s.(randLengthFluor_smooth)(f)=s.(randFluor_smooth)(f)*s.(randLength_smooth)(age);
              s.(randrealLengthFluor_smooth)(f)=s.(randFluor_smooth)(f)*s.(length_smooth)(age);
          end
   
          
          clear age;
        end
   end
   %------------------------------------------------------------------------
    schnitzcells(i) = s;
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES AND ADD RATES
%--------------------------------------------------------------------------  
for i = 1:length(schnitzcells)
  s = schnitzcells(i);  
  
  if length(s.(fluortime))==0 % e.g. weird long cells which divide very fast 
                              % have sometimes no fluor image
      continue  % (NW 2012-08)
  else
      %------------------------------------------------------------------------
      % LOOP OVER FLUOR MEASUREMENTS EXCEPT LAST
      %------------------------------------------------------------------------
      if length(s.(fluortime)) > 1
        for f = 1:length(s.(fluortime))-1
       %   age                 = find(s.time == s.(fluortime)(f)); % different to Daans version!
       %   age2                = find(s.time == s.(fluortime)(f+1)); % different to Daans version!

      %    s.(d_sum)(end+1)         =  s.(sum_all)(age2) - s.(sum_all)(age);
    %      s.(d_sum_dt)(end+1)      = (s.(sum_all)(age2) - s.(sum_all)(age)) /  (s.(fluortime)(f+1) - s.(fluortime)(f));
          s.(dFluor_mean_length_dt)(end+1)      = (s.(Fluor_mean_length)(f+1) - s.(Fluor_mean_length)(f)) /  (s.(fluortime)(f+1) - s.(fluortime)(f));
          s.(dFluor_mean_length_smooth_dt)(end+1)      = (s.(Fluor_mean_length_smooth)(f+1) - s.(Fluor_mean_length_smooth)(f)) /  (s.(fluortime)(f+1) - s.(fluortime)(f));
          if ADD_RND_DATA
                s.(drandLengthFluor_smooth_dt)(end+1)      = (s.(randLengthFluor_smooth)(f+1) - s.(randLengthFluor_smooth)(f)) /  (s.(fluortime)(f+1) - s.(fluortime)(f));
                s.(drandrealLengthFluor_smooth_dt)(end+1)      = (s.(randrealLengthFluor_smooth)(f+1) - s.(randrealLengthFluor_smooth)(f)) /  (s.(fluortime)(f+1) - s.(fluortime)(f));
          end

        end
      else
        f = 0;
      end
      %------------------------------------------------------------------------

      %------------------------------------------------------------------------
      % FOR LAST FLUOR MEASUREMENT CHECK WHETHER CHILDREN EXIST AND HAVE FLUOR
      % ADD DATA CURRENTLY ALSO FOR SMOOTHED CASE (might be changed)
      %------------------------------------------------------------------------
      clear sD sE;
      if (s.D ~= 0), sD = schnitzcells(s.D); end
      if (s.E ~= 0), sE = schnitzcells(s.E); end
      if (exist('sE') & exist('sD'))
        if (length(sD.(fluortime)) > 0 & length(sE.(fluortime)) > 0)
         if i==1 % enable e.g. if first schnitzcell has no fluor at all. Remove this schnitz in later analysis!
         %    
             s.(dFluor_mean_length_dt)         = 0;
             s.(dFluor_mean_length_smooth_dt)         = 0;  %BLUBB CAREFULL!!!
             s.(Fluor_mean_length_smooth)         = 0;  %BLUBB CAREFULL!!!
             s.(Fluor_mean_length)         = 0;  %BLUBB CAREFULL!!!
             if ADD_RND_DATA
                 s.(randLengthFluor_smooth) =0;  %BLUBB CAREFULL!!!
                 s.(drandLengthFluor_smooth_dt)=0; %BLUBB CAREFULL!!!
                 s.(drandrealLengthFluor_smooth_dt)=0; %BLUBB CAREFULL!!!
             end

         % s.(d_sum_dt)(end+1)      =  0;
         % s.(d_sum_dt_ph)(end+1)   = 0 ;
         % s.(d_sum_dt_vol)(end+1)  = 0;

        %  s.(phase_atfluor)(end+1)       =  0;
         else
       %   age                 = find(s.time == s.(fluortime)(end));
       %   age2                = find(sD.time == sD.(fluortime)(1));
      % i     
      % sD.(Fluor_mean_length)(1)
      %      sE.(Fluor_mean_length)(1)
      %      s.(Fluor_mean_length)(end) 
      %      sD.(fluortime)(1) 
      %      s.(fluortime)(end)


          s.(dFluor_mean_length_dt)(end+1)      = (sD.(Fluor_mean_length)(1) + sE.(Fluor_mean_length)(1) - s.(Fluor_mean_length)(end)) / (sD.(fluortime)(1) - s.(fluortime)(end));

          % *******CHANGE******************************!!!!!!!
          s.(dFluor_mean_length_smooth_dt)(end+1)      = (sD.(Fluor_mean_length_smooth)(1) + sE.(Fluor_mean_length_smooth)(1) - s.(Fluor_mean_length_smooth)(end)) / (sD.(fluortime)(1) - s.(fluortime)(end));
           % *******CHANGE******************************!!!!!!!

           if ADD_RND_DATA
              s.(drandLengthFluor_smooth_dt)(end+1)      = (sD.(randLengthFluor_smooth)(1) + sE.(randLengthFluor_smooth)(1) - s.(randLengthFluor_smooth)(end)) / (sD.(fluortime)(1) - s.(fluortime)(end));
              s.(drandrealLengthFluor_smooth_dt)(end+1)      = (sD.(randrealLengthFluor_smooth)(1) + sE.(randrealLengthFluor_smooth)(1) - s.(randrealLengthFluor_smooth)(end)) / (sD.(fluortime)(1) - s.(fluortime)(end));
           end

          end
        end
     end
      %------------------------------------------------------------------------

     schnitzcells(i) = s;
   end % exclude weird cells
end 

save(schnitzname,'schnitzcells');