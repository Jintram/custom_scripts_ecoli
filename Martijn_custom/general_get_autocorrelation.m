function [S, Sraw] = general_get_autocorrelation(y1)
% This function calculatates the autocorrelation for a vector with
% y-values. It can weigh certain parts of this vector more or less, and 
% will do so by weighing vector W. If w=[1,1,1,..,1], this function will
% give the same output as xcov(y1,'coeff').
%
% Note that S(tau=0) = S(1) b/c matlab starts indexing at 1

% First calculate distances to mean "mean removed y1".
mr_y1 = y1 - mean(y1);

% Give me 
% S(tau) = 1/((t_total-tau)) sum_0^(t_total-tau) [ (y1(tstar)y2(tstar+tau)) dtstar ]
% R(tau) = S(tau)/S(0)
% In discrete form (tau -> N*Dtau = NDtau) etc 
S = [];
for NDtau = 0:length(mr_y1)

    S_NDtau = [];
    
    % slide the window
    for NDtstar = 1:(length(mr_y1)-NDtau)
        
        % multiply two values
        S_NDtau(end+1) = mr_y1(NDtstar)*mr_y1(NDtstar+NDtau);
        
    end
    
    % calculate average and record
    S(NDtau+1) = sum(S_NDtau)/length(S_NDtau); % +1 b/c matlab starts indexing at 1
    Sraw(NDtau+1)= sum(S_NDtau);
    
end

end




