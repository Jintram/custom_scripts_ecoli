

%%
script20171031_Donachiedata

%%
figure;
%scatter(DonachieData(:,1), DonachieData(:,2),7^2,'filled');

scatter(DonachieDataPenicilin(:,2), DonachieDataPenicilin(:,1)./DonachieDataPenicilin(:,2),7^2,[.7 .7 .7],'filled');
hold on;
scatter(DonachieDataPenicilin(:,2), 1-DonachieDataPenicilin(:,1)./DonachieDataPenicilin(:,2),7^2,[.7 .7 .7],'filled');

ylim([0,1]);

title(['Donachie and Begg (1970)' 10 'Penicilin sensitive sites']);
xlabel('Cellular length');
ylabel('Relative position');
MW_makeplotlookbetter(20);

%%
figure;
%scatter(DonachieData(:,1), DonachieData(:,2),7^2,'filled');

scatter(DonachieDataDivisions(:,2), DonachieDataDivisions(:,1)./DonachieDataDivisions(:,2),7^2,[.7 .7 .7],'filled');
hold on;
scatter(DonachieDataDivisions(:,2), 1-DonachieDataDivisions(:,1)./DonachieDataDivisions(:,2),7^2,[.7 .7 .7],'filled');

ylim([0,1]);

title(['Donachie and Begg (1970)' 10 'Division sites penicilin recovery']);
xlabel('Cellular length');
ylabel('Relative position');
MW_makeplotlookbetter(20);