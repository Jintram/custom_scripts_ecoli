function Anorm = normalize(A)
    % Rescales matrix such that values go from 0-1.   
    % - MW 2014/05
    
    Anorm = (A - min(A(:))) ./ (max(A(:)) - min(A(:)));