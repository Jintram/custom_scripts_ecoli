% Principal Component Analysis

% **** ADJUST *****
schnitzUseName='schnitzcells_malt20120508pos2_full';  %schnitzcells selection
%field1='dY5_sum_dt_subtr';  
%field2='dC5_sum_dt_subtr';
%field3='muP15_fitNew_subtr'; 
%timefield='time_atdY'; % timefield of field1,2 (rates) maybe shorter than of field3 (growth rate).
%                         But t has to start with the same time (no shift)
%                         and cannot be longer
%field1='dY5_sum_dt_s_cycCor';  
%field2='dC5_sum_dt_s_cycCor';
%field3='muP15_fitNew_cycCor'; 
%timefield='time_atdY';
field1='dY5_cycCor';  
field2='dC5_cycCor';
field3='muP15_fitNew_cycCor'; 
timefield='dY5_time';

PLOTPARETO=0;
PLOTBIPLOT=1;
% *****************

% get schitzcells structure (above: only name)
eval(['schnitzUse=' schnitzUseName ';']);

% extract data
mymatrix=zeros(0,3);

% loop over all schnitzes
for i=1:length(schnitzUse)
    s=schnitzUse(i);
    % use only schnitzes with useForPlot=1
    if s.useForPlot==1
        % check if fields actually contain data (length>0)
        % Note: e.g. prod.rates and mu do not always have same array length
        % (rates can be 1 shorter for late data)
        minLength=min([length(s.(field1)),length(s.(field2)),length(s.(field3))]);
        if length(s.(timefield))~=minLength
            error('Something wrong with timefield and length of PCA fields')
        end
        if minLength>0
            % check if NaN values exist
            Nan1=sum(isnan(s.(field1)(1:minLength)));
            Nan2=sum(isnan(s.(field2)(1:minLength)));
            Nan3=sum(isnan(s.(field3)(1:minLength)));
            if (Nan1+Nan2+Nan3)==0
                % add values to mymatrix
                addmatrix=[s.(field1)(1:minLength)',s.(field2)(1:minLength)',s.(field3)(1:minLength)'];
                mymatrix=[mymatrix;addmatrix];
            end
        end
    end
end

            
%%
% (1) *** normalize mean and std: ***
usematrix=zscore(mymatrix);
% (1end)
% (2)
% *** normalize mean : ***
%usematrix=mymatrix;
%for i=1:3
%    usematrix(:,i)=usematrix(:,i)./mean(usematrix(:,i));
%end
% (2end)

[coefs,scores,variances,t2] = princomp(usematrix);
disp(['Analyzed Princ Comp of ' schnitzUseName]);

coefs
variances=variances'
percent_explained = 100*variances/sum(variances)'
sumpercent_explained = cumsum(percent_explained)'
if PLOTPARETO
    figure
    pareto(percent_explained)
    xlabel('Principal Component')
    ylabel('Variance Explained (%)')
    title(['Eigenvalues of PCA for ' schnitzUseName],'Interpreter','none');
end
%%
%1:2 instead of 1:3 : 2 dimensional projection on 1st&2nd PCA component
if PLOTBIPLOT
    figure
    set(gcf, 'defaulttextinterpreter','none')
    clf
    biplot(coefs(:,1:3),'scores',scores(:,1:3),'varlabels', {field1 field2 field3} )
    title(['biplot of PCA for ' schnitzUseName]);
end
%%
%figure
%plot3(scores(:,1),scores(:,2),scores(:,3),'.g')
%xlabel('Component 1')
%ylabel('Component 2')
%zlabel('Component 3')
%%
%figure
%plot3(-usematrix(:,1),-usematrix(:,2),-usematrix(:,3),'.g')
%xlabel(field1)
%ylabel(field2)
%zlabel(field3)
