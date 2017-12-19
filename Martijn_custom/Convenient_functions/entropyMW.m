function H=entropyMW(pdf)
% Calculates the entropy of a probability distribution function.
% H = -sum(pdfNormalized.*log2(pdfNormalized))
%
% Normalizes the function for you, and also can handle 2d functions.
% Entries with a value smaller than EPSILON will be ignored (EPSILON is a
% parameter set in the function).
EPSILON=10e-14;

% normalize
pdfNormalized = pdf./sum(pdf(pdf>EPSILON));

% calculate entropy (take into account support only)
H = -sum(pdfNormalized(pdfNormalized>EPSILON).*log2(pdfNormalized(pdfNormalized>EPSILON)));
    % https://en.wikipedia.org/wiki/Differential_entropy
    % https://en.wikipedia.org/wiki/Support_(mathematics)

%{
if size(pdfNormalized,1)>1 & size(pdfNormalized,2)>1
    % 2d
    H = -sum(sum(pdfNormalized.*log2(pdfNormalized)))
else
    % 1d
    H = -sum(pdfNormalized.*log2(pdfNormalized))
end
%}

end