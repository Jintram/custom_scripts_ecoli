

% Joris' formula

m=[1:1000];
N=2000;
fraction=.5;

figure, plot(m, m./N + (1./m - 1./N) .* fraction,'o');

% Note that at m = 1, the fraction is .2. But the line for high m is simply
% dominated by m/N behavior.. 

% Unfortunately this all didn't work..