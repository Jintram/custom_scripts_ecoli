

%% 

% Note that below values come from 
% mkate2
% venus 
% cerulean
% TurboGFP
name_fp   = {'mKate'  , 'mVenus', 'mCerulean', 'GFP'};
lambda_ex = [588      , 515     , 433        , 482  ];
lambda_em = [633      , 528     , 475        , 502  ];

myColors = [255,66, 0; 86,255, 0; 0,192, 255;0,255, 123]./255; % emission converted with https://academo.org/demos/wavelength-to-colour-relationship/

f(1)=figure(1); clf; hold on;
scatter(lambda_ex(1:3), lambda_em(1:3),15^2,myColors(1:3,:),'filled','MarkerEdgeColor','k');
scatter(lambda_ex(4), lambda_em(4),15^2,myColors(4,:),'filled','s','MarkerEdgeColor',[.7 .7 .7]);
xlim([350, 750]);
ylim([400, 750]);

for idx=1:4
    text(lambda_ex(idx)+15,lambda_em(idx),name_fp{idx});
end

xlabel('Excitation wavelength');
ylabel('Emission wavelength');

MW_makeplotlookbetter(10*2);%,optionalParameters);

%%
% Note that our filters are: Chroma, 41017, 49008, 49001 and 49003
% (resp. GFP, mCherry, CFP and YFP)
% Source of data:
% https://www.chroma.com/products/sets/41017-endow-gfp-egfp-bandpass#tabs-0-main-1
% and other sets

% Import data (see below for names)
mcherry_ex = importdata('D:\Local_Software\Martijn_extensions\Martijn_custom\misc\data\mcherry_ex.txt');
mcherry_em = importdata('D:\Local_Software\Martijn_extensions\Martijn_custom\misc\data\mcherry_em.txt');

eyfp_ex= importdata('D:\Local_Software\Martijn_extensions\Martijn_custom\misc\data\eyfp_ex.txt');
eyfp_em= importdata('D:\Local_Software\Martijn_extensions\Martijn_custom\misc\data\eyfp_em.txt');

ecfp_ex= importdata('D:\Local_Software\Martijn_extensions\Martijn_custom\misc\data\ecfp_ex.txt');
ecfp_em= importdata('D:\Local_Software\Martijn_extensions\Martijn_custom\misc\data\ecfp_em.txt');

engfp_ex= importdata('D:\Local_Software\Martijn_extensions\Martijn_custom\misc\data\engfp_ex.txt');
engfp_em= importdata('D:\Local_Software\Martijn_extensions\Martijn_custom\misc\data\engfp_em.txt');

%% Plot excitation in filters
f(2)=figure(2); clf; hold on;

plot(mcherry_ex(:,1), mcherry_ex(:,2),'Color',myColors(1,:),'LineWidth',2);
plot(eyfp_ex(:,1), eyfp_ex(:,2),'Color',myColors(2,:),'LineWidth',2);
plot(ecfp_ex(:,1), ecfp_ex(:,2),'Color',myColors(3,:),'LineWidth',2);
plot(engfp_ex(:,1), engfp_ex(:,2),'--','Color',myColors(4,:),'LineWidth',2);

% add excitation values proteins
plot(lambda_ex(1),0,'^','Color','k','MarkerSize',15,'MarkerFaceColor',myColors(1,:));
plot(lambda_ex(2),0,'^','Color','k','MarkerSize',15,'MarkerFaceColor',myColors(2,:));
plot(lambda_ex(3),0,'^','Color','k','MarkerSize',15,'MarkerFaceColor',myColors(3,:));
plot(lambda_ex(4),0,'^','Color',[.7 .7 .7],'MarkerSize',10,'MarkerFaceColor',myColors(4,:));

xlabel('Excitation wavelength');
ylabel('Transmission');
xlim([350, 750]);

MW_makeplotlookbetter(10*2);%,optionalParameters);

%% Plot emission in filters
f(3)=figure(3); clf; hold on;

plot(mcherry_em(:,2), mcherry_em(:,1),'Color',myColors(1,:),'LineWidth',2);
plot(eyfp_em(:,2), eyfp_em(:,1),'Color',myColors(2,:),'LineWidth',2);
plot(ecfp_em(:,2), ecfp_em(:,1),'Color',myColors(3,:),'LineWidth',2);
plot(engfp_em(:,2), engfp_em(:,1),'--','Color',myColors(4,:),'LineWidth',2);

% add emission values proteins
plot(0,lambda_em(1),'>','Color','k','MarkerSize',15,'MarkerFaceColor',myColors(1,:));
plot(0,lambda_em(2),'>','Color','k','MarkerSize',15,'MarkerFaceColor',myColors(2,:));
plot(0,lambda_em(3),'>','Color','k','MarkerSize',15,'MarkerFaceColor',myColors(3,:));
plot(0,lambda_em(4),'>','Color',[.7 .7 .7],'MarkerSize',10,'MarkerFaceColor',myColors(4,:));


xlabel('Emission wavelength');
ylabel('Transmission');
ylim([400, 750]);

MW_makeplotlookbetter(10*2);%,optionalParameters);

%%

OUTPUTFOLDER = '\\storage01\data\AMOLF\users\wehrens\Latex3\Thesis\Chapter2_Methods\Figures\MatlabExport\';

subLabel={'a','b','c'};

for fIdx=1:3

    currentf=f(fIdx);

    figure(currentf);
    SIZE=[6.80,6.80]; OFFSET = [0,0];
    set(currentf,'Units','centimeters','Position',[OFFSET SIZE]*2)
    MW_makeplotlookbetter(10*2);%,optionalParameters);
    set(currentf,'RendererMode','manual','Renderer','Painters');

    filename= ['fluorlabels_' subLabel{fIdx}];
    saveas(currentf,[OUTPUTFOLDER filename '.svg']); saveas(currentf,[OUTPUTFOLDER filename '.tif']); saveas(currentf,[OUTPUTFOLDER filename '.fig']);
    
end
