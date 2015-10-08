
%% Generate some testdata
testData = rand(1,100000);

%% make a histogram
[elements, centers] = hist(testData); % create histogram
dx = centers(2)-centers(1); % for normalization later
area = sum(elements)*dx; % for normalization later

% make a KDE estimate
[f,xi]=ksdensity(testData)

%% plot it
figure(1); clf; hold on;
plot( centers, elements./area,'o','LineWidth',2) % histogram
plot(xi,f,'--k','LineWidth',2) % KDE estimate

% appareance
ylim([0,2]);
xlabel('x')
ylabel('pdf(x)')

% fontsize of plot
FONTSIZE=20
set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal')
set(gca,'FontSize',FONTSIZE)