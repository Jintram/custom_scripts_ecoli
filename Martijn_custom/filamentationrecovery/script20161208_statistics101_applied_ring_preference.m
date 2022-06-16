


% Simple statistics on rings


% probability mass function:
% f(k,n,p) = (n./k) .* p.^k (1-p).^(n-k)
% https://en.wikipedia.org/wiki/Binomial_distribution

pmf = @(k,n,p) nchoosek(n,k) .* p.^k .* (1-p).^(n-k);

pmf_arrayfun = @(k,n,p) arrayfun(@(k) pmf(k,n,p), k);


n=51;% 288;
k=[1:n];
p=0.66;
pmf_values = pmf_arrayfun(k,n,p);

figure; plot(k,pmf_values);