function R = general_get_crosscorrelation(X)
% This function calculatates the crosscorrelation between two vectors with
% y-values. It can weigh certain parts of this vector more or less, and 
% will do so by weighing vector W. If w=[1,1,1,..,1], this function will
% give the same output as xcorr.

% Give me 
% S(tau) = integral_0^(N-Ntau) [ (y1(tstar)y2(tstar+tau)) dtstar ]
% R(tau) = S(tau)/S(0)
for tau = 1:length(X)

    for tstar = 1:(length(X)-tau)
        
        
        
    end
    
end

end