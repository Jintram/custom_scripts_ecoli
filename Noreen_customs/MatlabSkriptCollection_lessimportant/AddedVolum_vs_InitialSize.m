Ld=[];
Lb=[];
avmu=[];
tau=[];

ss=schnitzcells_rm;
for i=1:length(ss)
    if ss(i).useForPlot==1 & ss(i).completeCycle==1
        Ld=[Ld ss(i).length_fitNew(end)];
        Lb=[Lb ss(i).length_fitNew(1)];
        avmu=[avmu ss(i).av_mu_fitNew ];
        tau=[tau ss(i).interDivTime];
    end
end

DeltaL=Ld-Lb;