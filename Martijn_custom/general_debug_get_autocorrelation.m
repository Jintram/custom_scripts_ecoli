
%% Check weights

mysums = [];

for NDtau = 0:length(myY)
    mysums(end+1) = sum(w(:,NDtau+1));
end

mysums

%% Check whether dt_i == dt_j
dts = [];
for idx = 1:length(myti)-1
    dts(end+1) = myti(idx+1)-myti(idx);
end

dts