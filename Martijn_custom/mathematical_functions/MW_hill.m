function y = MW_hill(myParameters, L)

n = myParameters(1);
K = myParameters(2);
ymax = myParameters(3);

y = ymax * L.^n ./ (L.^n + K);

end