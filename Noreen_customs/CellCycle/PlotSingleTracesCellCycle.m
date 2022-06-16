%schnitznumrange=[49 61 102 142 167 182 210 255 304 312 32  354 369 380 388 422 433 446 463];
%schnitznumrange=allstepsVec(:,1)';
schnitznumrange=[200:220];
myfield='dY5';   % alternatively dC5, dY5_sum_dt etc
myphase='phase2_at_dY5_time';
PLOTCONC=1;
concfield='Y5_mean';
concphase='phase2_atY';
schnitzUse=s_rm_fitTime;
figure
for schn=schnitznumrange
    if schnitzUse(schn).useForPlot==1 & schnitzUse(schn).completeCycle==1
        yy=schnitzUse(schn).(myfield);
        phph=schnitzUse(schn).(myphase);
        if ~PLOTCONC
            figure
            title([num2str(schn)]);
            hold on
            plot(phph,yy,'LineWidth',2)
        else
            yyconc=schnitzUse(schn).(concfield);
            phphconc=schnitzUse(schn).(concphase);
            %figure
            hold on
            subplot(3,1,1)
            plot(phph,yy)
            title([num2str(schn)]);
            xlabel('rate')
            subplot(3,1,2)
            plot(phphconc,yyconc,'r')
            xlabel('conc')
            subplot(3,1,3)
            schnitzlength=schnitzUse(schn).length_fitNew;
            schnitzphase=schnitzUse(schn).phase;
            plot(schnitzphase,schnitzlength,'k')
            xlabel('length')
        end
    end
end




%%
schnitznumrange=[400:420 ];
myfield='dR5';   % alternatively dC5, dY5_sum_dt etc
myphase='phase2_at_dR5_time';
schnitzUse=s_rm_fitTime;

for schn=schnitznumrange
    if schnitzUse(schn).useForPlot==1 & schnitzUse(schn).completeCycle==1
        yy=schnitzUse(schn).(myfield);
        phph=schnitzUse(schn).(myphase);
        figure
        title([num2str(schn)]);
        hold on
        plot(phph,yy/mean(yy),'r')
        plot(phph,schnitzUse(schn).dG5/mean(schnitzUse(schn).dG5),'g')
    end
end