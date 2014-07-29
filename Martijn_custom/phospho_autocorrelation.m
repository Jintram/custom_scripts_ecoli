
% So we want an autocorrelation, e.g. of the following data:
myT = current_branchData(bac_nr).frame_nrs;
myY = current_branchData(bac_nr).muP11_all;

% get manual ac
[ac,ac_raw]=general_get_autocorrelation(myY);

% get matlab ac
[ac_ml,lags]=xcov(myY,'coeff');

figure(1);
clf; hold on;
plot(ac/ac(1),'-b');
plot(ac_raw/ac_raw(1),'--b','LineWidth',2);
plot(lags,ac_ml,'-r');