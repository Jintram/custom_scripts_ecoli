
branch_nr = 200;

% So we want an autocorrelation, e.g. of the following data:
myT = current_branchData(branch_nr).frame_nrs;
myY = current_branchData(branch_nr).muP11_all;
myLambda = current_branchData(branch_nr).count; % number of time datapoint is used for in branch structure

% get manual ac without correction for branches
[ac,ac_raw]=general_get_autocorrelation(myY,0);
ac_w = general_get_autocorrelation(myY,myLambda);

% get matlab ac
[ac_ml,lags]=xcov(myY,'coeff');

figure(1);
clf; hold on;
plot(ac/ac(1),'-b');
plot(ac_raw/ac_raw(1),'--b','LineWidth',2);
plot(ac_w/ac_w(1),'--g','LineWidth',2);

plot(lags,ac_ml,'-r');
