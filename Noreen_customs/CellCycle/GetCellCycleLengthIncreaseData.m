%% Total length increase
ph=[];
totlength=[];
timedata=[];
schn=[];
frames=[];
frames_unique=[];
highlightschnitzes=[111 200 444];
highlightcolors=[0.6 0 0; 0.6 0 0; 0.6 0 0]; % must have the same #columns as length of highlightschnitzes
RESTRICTTOFLUOIMAGES=0;  % only take data which appears at fluorphasefield
fluoframesall='Y_frames_all'; % concentration/fluo image time points. late cells may have one fluorate datapoint less
                           % when no fluo image exists, value is 'nan'

myschnitzcells=s_rm_fitTime;

for i=1:length(myschnitzcells)
    s=myschnitzcells(i);
    if s.useForPlot==1
        if s.completeCycle==1 & ~isempty(s.time) & length(s.phase)==length(s.length_fitNew) & length(s.phase)==length(s.time) 
            if ~RESTRICTTOFLUOIMAGES %take all data
                ph=[ph; s.phase'];
                totlength=[totlength; s.length_fitNew']; %smoothing?            
                %totlength=[totlength; s.length_fitNew'/mean(s.length_fitNew)];
                             % normalize by average
                %totlength=[totlength; s.length_fitNew'/(s.length_fitNew(1))];
                            % normalize by initial value
                reltime=s.time-s.time(1);
                timedata=[timedata; reltime'];
                schnrep=zeros(1,length(s.phase))+i;
                schn=[schn; schnrep'];
                framesrep=zeros(1,length(s.phase))+length(s.frames);
                frames=[frames; framesrep'];
                frames_unique=[frames_unique; length(s.frames)];
            else %take only data at fluotimefields
                fluoframeidx=find(~isnan(s.(fluoframesall)));
                ph=[ph; s.phase(fluoframeidx)'];
                totlength=[totlength; s.length_fitNew(fluoframeidx)']; %smoothing?            
                %totlength=[totlength; s.length_fitNew(fluoframeidx)'/mean(s.length_fitNew(fluoframeidx))];
                            % normalize by average
                %totlength=[totlength; s.length_fitNew(fluoframeidx)'/(s.length_fitNew(fluoframeidx(1)))];
                            % normalize by initial value
                reltime=s.time(fluoframeidx)-s.time(1);
                timedata=[timedata; reltime'];
                schnrep=zeros(1,length(s.phase(fluoframeidx)))+i;
                schn=[schn; schnrep'];
                framesrep=zeros(1,length(s.phase(fluoframeidx)))+length(s.frames);
                frames=[frames; framesrep'];
                frames_unique=[frames_unique; length(s.frames)];
            end   
        end
    end
end

% normalize and save unnormalized vector
totlengthabs=totlength;
totlength=totlength/mean(totlength);

% binning LENGTH
%idx=find(frames<61 & frames>54);
idx=find(frames<665 & frames>0);
%plot(ph(idx),totlength(idx),'.r')
phsub=ph(idx);
totlengthsub=totlength(idx);


% do binning
numbins=10;
binborders=[0:1/numbins:0.98];
binborders=[binborders,1];
binnedtotlength=[];
meanbin=[0.5/numbins:1/numbins:0.98];

for i=1:numbins
    idx2=find(phsub>=binborders(i) & phsub<binborders(i+1));
    binnedtotlength=[binnedtotlength,mean(totlengthsub(idx2))];
end

% fit exponential to bins
[expo, L0]=DJK_ExponentialFit(meanbin,binnedtotlength);
phasecont=0:0.01:1;
lengthcont=L0*2.^(expo*phasecont);

%binned length increase
figure(2)
clf
plot(meanbin,binnedtotlength,'.-','MarkerSize',15,'Color','r')
hold on
plot(phasecont, lengthcont,'k')
grid on


%highlightschnitzes=[444 555 666];
figure(1)
clf
hold on
% ***** plot dots
plot(ph(idx),totlength(idx),'.','Color','r')
%plot lines
%for j=1:length(unique(schn))
%    schnunique=unique(schn);
%    ss=schnunique(j);
%    idxschn=find(schn==ss);
%    phschn=ph(idxschn);
%    totlengthschn=totlength(idxschn);
%    plot(phschn,totlengthschn,'-','Color','r')
%end
% *****
plot(meanbin,binnedtotlength,'.-k','MarkerSize',15)
if ~isempty(highlightschnitzes)
    if size(highlightcolors,1)==length(highlightschnitzes);
        for j=1:length(highlightschnitzes)
            ss=highlightschnitzes(j);
            idxschn=find(schn==ss);
            phschn=ph(idxschn);
            totlengthschn=totlength(idxschn);
            plot(phschn,totlengthschn,'.-','Color',highlightcolors(j,:))
        end
    else
        disp('wrong # of highlighted lineage colours')
    end
end
    


