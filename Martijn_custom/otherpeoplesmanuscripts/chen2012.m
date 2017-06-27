

% extracted from chen2012 using http://arohatgi.info/WebPlotDigitizer/app/
data30Slocations = 14:45
data30S =...
    [[0, 0.10224438902743153];...
    [1, 0.35411471321695764];...
    [2, 0.4438902743142144];...
    [3, 0.4139650872817955];...
    [4, 0.3690773067331671];...
    [5, 0.32668329177057365];...
    [6, 0.29675810473815467];...
    [7, 0.25436408977556113];...
    [8, 0.23441396508728188];...
    [9, 0.2443890274314215];...
    [10, 0.34413965087281806];...
    [11, 0.40149625935162103];...
    [12, 0.28428927680798005];...
    [13, 0.19451371571072332];...
    [14, 0.13965087281795516];...
    [15, 0.13216957605985036];...
    [16, 0.11471321695760597];...
    [17, 0.11720698254364094];...
    [18, 0.1421446384039901];...
    [19, 0.2493765586034913];...
    [20, 0.3640897755610973];...
    [21, 0.2768079800498753];...
    [22, 0.2618453865336659];...
    [23, 0.28428927680798005];...
    [24, 0.41645885286783046];...
    [25, 0.6708229426433915];...
    [26, 0.7905236907730674];...
    [27, 0.5561097256857854];...
    [28, 0.2793017456359103];...
    [29, 0.17206982543640895];...
    [30, 0.11720698254364094];...
    [31, 0.08478802992518714]];

%{
    data30Slocations = 20:45
    [[0, 0.024999999999999887];...
    [1, 0.11666666666666656];...
    [2, 0.2666666666666666];...
    [3, 0.19999999999999993];...
    [4, 0.09999999999999996];...
    [5, 0.04999999999999998];...
    [6, 0.029166666666666535];...
    [7, 0.024999999999999887];...
    [8, 0.029166666666666535];...
    [9, 0.03749999999999983];...
    [10, 0.05416666666666663];...
    [11, 0.07499999999999987];...
    [12, 0.09166666666666667];...
    [13, 0.1375];...
    [14, 0.39999999999999997];...
    [15, 0.5208333333333334];...
    [16, 0.3083333333333333];...
    [17, 0.16041666666666657];...
    [18, 0.10833333333333325];...
    [19, 0.04999999999999998];...
    [20, 0.04583333333333334];...
    [21, 0.05833333333333328];...
    [22, 0.09999999999999996];...
    [23, 0.06249999999999993];...
    [24, 0.04583333333333334];...
    [25, 0.05833333333333328]];
%}
%{
[[0, 0.020833333333333242];...
[1, -0.004166666666666851];...
[2, 0.008333333333333297];...
[3, 0.26249999999999996];...
[4, 0.19999999999999993];...
[5, 0.1041666666666666];...
[6, 0.03749999999999983];...
[7, 0.024999999999999887];...
[8, 0.024999999999999887];...
[9, 0.03333333333333319];...
[10, 0.03749999999999983];...
[11, 0.05416666666666663];...
[12, 0.07083333333333322];...
[13, 0.09583333333333333];...
[14, 0.13333333333333336];...
[15, 0.3958333333333333];...
[16, 0.5041666666666667];...
[17, 0.31249999999999994];...
[18, 0.15833333333333324];...
[19, 0.1041666666666666];...
[20, 0.04999999999999998];...
[21, 0.04999999999999998];...
[22, 0.05833333333333328];...
[23, 0.09583333333333333];...
[24, 0.04583333333333334];...
[25, 0.05833333333333328]];
    %}
    
data30S = data30S-min(data30S(:));

figure(1); clf; hold on;
bar(data30Slocations,data30S(:,2)');

xlabel('Fraction');
ylabel('Abundance');
MW_makeplotlookbetter(20);

%%

S2  = [ 0 0 0.13 0 0.11 0.17 0.24 0.71 0.76 0.55 0.35 0.35 0.39 0.34 0.38 0.48 0.41 0.22 0.38 0.49 0.76 0.67 0.55 0.44 0.56 0.56 0.71 0.75 0.71 0.61 0.70 0.74 ];
S3  = [ 0.34 0.36 0.22 0.15 0.10 0.17 0.24 0.72 0.76 0.55 0.35 0.34 0.38 0.33 0.36 0.46 0.40 0.21 0.37 0.48 0.76 0.67 0.54 0.44 0.55 0.56 0.70 0.74 0.70 0.60 0.70 0.74];
S4  = [ 0.29 0.36 0.44 0.43 0.42 0.45 0.39 0.76 0.79 0.58 0.38 0.38 0.42 0.37 0.40 0.49 0.42 0.22 0.37 0.48 0.76 0.67 0.55 0.43 0.56 0.56 0.71 0.74 0.71 0.60 0.70 0.73];
S14 = [ 0.14 0.12 0 0.12 0.12 0.17 0.26 0.74 0.78 0.56 0.37 0.35 0.41 0.35 0.37 0.47 0.42 0.22 0.38 0.48 0.76 0.67 0.54 0.43 0.56 0.56 0.71 0.74 0.70 0.60 0.69 0.73];

%{
S2=...
    [0.24 0.71 0.76 0.55 0.35 0.35 0.39 0.34 0.38 0.48 0.41 0.22 0.38 0.49 0.76 0.67 0.55 0.44 0.56 0.56 0.71 0.75 0.71 0.61 0.70];
S3=...
    [0.24 0.72 0.76 0.55 0.35 0.34 0.38 0.33 0.36 0.46 0.40 0.21 0.37 0.48 0.76 0.67 0.54 0.44 0.55 0.56 0.70 0.74 0.70 0.60 0.70 0.74];
S4 = ...
    [0.39 0.76 0.79 0.58 0.38 0.38 0.42 0.37 0.40 0.49 0.42 0.22 0.37 0.48 0.76 0.67 0.55 0.43 0.56 0.56 0.71 0.74 0.71 0.60 0.70 0.73];
S14 = ...
    [0.26 0.74 0.78 0.56 0.37 0.35 0.41 0.35 0.37 0.47 0.42 0.22 0.38 0.48 0.76 0.67 0.54 0.43 0.56 0.56 0.71 0.74 0.70 0.60 0.69 0.73];
%}

figure(2); clf; hold on;

plot(S2,'o-','LineWidth',2);
plot(S3,'o-','LineWidth',2);
plot(S4,'o-','LineWidth',2);
plot(S14,'o-','LineWidth',2);

xlabel('Fraction');
ylabel('Protein abundance');
legend({'S2','S3','S4','S14'});
MW_makeplotlookbetter(20);

%% Normalized plot

totalpdf = sum(data30S(:,2)');

normS3=S3     .* data30S(:,2)'./totalpdf;
normS4=S4     .* data30S(:,2)'./totalpdf;
normS14=S14   .* data30S(:,2)'./totalpdf;

figure(3); clf; hold on; title('Normalized');
%plot(S2.*data30S(:,2)');
plot(data30Slocations,normS3,'o-','LineWidth',2);
plot(data30Slocations,normS4,'o-','LineWidth',2);
plot(data30Slocations,normS14,'o-','LineWidth',2);

maxS14=max(normS14(~isinf(normS14)))
ylim([0,maxS14*1.1]);

xlabel('Fraction');
ylabel('Protein abundance');
legend({'S3 (late)','S4 (early)','S14'});
MW_makeplotlookbetter(20);


%%

