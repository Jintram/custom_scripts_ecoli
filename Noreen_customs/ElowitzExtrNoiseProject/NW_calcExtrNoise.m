function [noise_extr,noise_intr]=NW_calcExtrNoise(schnitzcells,field1,field2)
% extrinsic and intrinsic noise according to formulas in Elowitz2002 are
% calculated.
%
% Formulas used: noise_int^2=<(c-y)^2>/(2<c><y>)
%                noise_ext^2=(<cy>-<c><y>)/(<c><y>)=corr(c,y)
%                noise_tot^2=(<c^2+y^2>-2<c><y>)/(2<c><y>)=noise_int^2+noise_ext^2
%
%
% field1, field2: schnitzcells fields for which noise is calculated.
%                 They have to be one of these pairs (!!!):
%                 (Y6_mean, C6_mean)
%                 (YFPrate, CFPrate)
%                 maybe: (av_Y6_mean, av_C6_mean)
%      No check if wrong combination (e.g. conc&rate) was used!!!
%
% uses only schnitzes with useForPlot==1. No weighing performed.

% there might be a former version 'calcExtrNoise' flying around

% *************** CHECK CODE UND NORM!!!

Y=[];
C=[];
YC=[];

for i=1:length(schnitzcells)
    if schnitzcells(i).useForPlot==1
        if length(schnitzcells(i).(field1))>0
            Y=[Y, schnitzcells(i).(field1)];
            C=[C, schnitzcells(i).(field2)];
            YC=[YC, schnitzcells(i).(field1).*schnitzcells(i).(field2)];%unnecessary
        end
    end
end
% remove nan values here if necessary... idx=~isnan(...)

% Signals have to be normalized so that they have same mean -> use mean=1
% (otherwise, e.g. (yfp-cfp) cannot stand for intrinsic component but is
% just fluor protein specific)
%norm=1;
%if norm==1
Y=Y/mean(Y);
C=C/mean(C);
YC=Y.*C;
%end
    
YY=mean(Y);
CC=mean(C);
YYCC=mean(YC);

noise_ext=sqrt((YYCC-YY*CC)/(YY*CC));

noise_int=sqrt(mean((C-Y).^2)/(2*CC*YY));

noise_total=sqrt((mean(C.^2+Y.^2)-2*CC*YY)/2*CC*YY);

% some outout
disp(['Field1: ' field1   '.   noise_ext=' num2str(noise_ext) '   noise_int=' num2str(noise_int)...
    '  noise_total=' num2str(noise_total)])
