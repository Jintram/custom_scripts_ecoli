% chooses schnitzes which start AND divide within fitTimeComplete
% output: s_rm_fitTime_birthdiv

%fitTimeComplete=[200 650];
fitTimeComplete=[100 700];


s_rm_fitTime_birthdiv=s_rm;

for i=1:length(s_rm)
    birth=s_rm_fitTime_birthdiv(i).birthTime;
    divi=s_rm_fitTime_birthdiv(i).divTime;
    % if no complete cell cycle
    if isnan(birth) | isnan(divi)
       s_rm_fitTime_birthdiv(i).useForPlot=0;
    % if cell cycle not within fitTime
    else
        mintime=fitTimeComplete(1);
        maxtime=fitTimeComplete(2);
        if birth<mintime | birth>maxtime |  divi<mintime | divi>maxtime
            s_rm_fitTime_birthdiv(i).useForPlot=0;
        end
    end
end



%%

    