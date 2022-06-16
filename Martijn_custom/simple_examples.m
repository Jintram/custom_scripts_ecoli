


%% Michaelis Menten curve
% Simple example

C = [0:.1:100];
Vmax = 3;
halfmax = 5;
V = (C.*Vmax)./(C+halfmax);

figure(1), plot(C,V);