% adds uniformly distributed random fields with length #frames and length
% #fluoframes to each schnitz.
%neames: rand_len, rand_fluor
%
function [schnitzcells] = NW_addToSchnitzes_random_Fields(p) 

schnitzname = [p.tracksDir,p.movieName,'-Schnitz.mat'];
load(schnitzname);

for i=1:length(schnitzcells)
    schnitzcells(i).rand_len=[];
    schnitzcells(i).rand_fluor=[];
    
    s=schnitzcells(i);
  s.rand_len=rand(1,length(s.time));
  s.rand_fluor=rand(1,length(s.Y_time));  %specific to YFP!
  
  schnitzcells(i)=s;
end 

save(schnitzname,'schnitzcells');