

myFileListing=dir('\\ASPARAGINE2\BigData\MartijnW\2014-10-10\pos14\')

dates=[]
deltas=[]
counter=0
for theFile=myFileListing'
    if length(theFile.name)>9
    if strcmp(theFile.name(1:9),'pos14-p-3')
        counter=counter+1;
        dates(end+1)=datenum(theFile.date);
        if counter>1
            deltas(end+1)=dates(end)-dates(end-1);
        end
    end
    end
end    
    
[counts, centers] = hist(deltas,100);
rescaling=24*60;
centersrescaled=centers.*rescaling;
figure(1), plot(centersrescaled, counts);
xlim([0,0.05].*rescaling);
%ylim([0,10]);
xlabel('delay between pictures (minutes)');
ylabel('count (#)');

