%%

FIELDNAME1 = 'G_time'
FIELDNAME2 = 'G6_mean_cycCor'
FIELDNAME3 = 'schnitzNrs'

count1=[]
count2=[]
count3=[]
for i = 1:numel(branches)
    count1 = [count1 numel(branches(i).(FIELDNAME1))]
    count2 = [count2 numel(branches(i).(FIELDNAME2))]
    count3 = [count3 numel(branches(i).(FIELDNAME3))]
end

%%

mintimes=[]; maxtimes=[]
for i = 1:numel(branches)
    mintimes = [mintimes min( branches(i).(FIELDNAME1) )]
    maxtimes = [maxtimes max( branches(i).(FIELDNAME1) )]
end

%% 
any(count1<64)