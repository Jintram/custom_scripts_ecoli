function [S, NDtau] = general_get_autocorrelation(y1,lambda,settings)
% This function calculatates the autocorrelation for a vector with
% y-values. It can weigh certain parts of this vector more or less, and 
% will do so by weighing vector W. If w=ones(N_ttotal), this function will
% give the same output as xcov(y1,'coeff').
%
% Note that S(tau=0) = S(1) b/c matlab starts indexing at 1
%
% If w=0, then w will be set to ones(N)
%
% input
% y1    
%    vector to calculate autocorr for
% lambda 
%    weighing factor, lambda=0 no weighing
% settings
%    additional options, settings=struct (ie empty) no settings
%

% Some general stuff
% ===
% Total nr points
N=length(y1); 
% Importantly, use distances to mean ("mean removed y1").
meanremoved_y1 = y1 - mean(y1); % 

if isfield(settings,'maxtau')
    maxtau = settings.maxtau;
else
    maxtau=length(meanremoved_y1);
end  

% Actual correlation
% ===
% Give me 
% S(tau) = 1/((t_total-tau)) sum_0^(t_total-tau) [ (y1(tstar)y2(tstar+tau)) dtstar ]
% R(tau) = S(tau)/S(0)
% In discrete form (tau -> N*Dtau = NDtau) etc 
S = []; 
for NDtau = 0:maxtau

    S_NDtau = [];
    %ws = []; % TEST
    
    % slide the window
    for NDtstar = 1:(length(meanremoved_y1)-NDtau)
        
        % multiply two values
        S_NDtau(end+1) = meanremoved_y1(NDtstar)*meanremoved_y1(NDtstar+NDtau);
        
        % multiply two values TESTING
        %{
        if lambda==0, lambdat1=1; else lambdat1=lambda(NDtstar); end
        if lambda==0, lambdat2=1; else lambdat2=lambda(NDtstar+NDtau); end
        S_NDtau(end+1) = (mr_y1(NDtstar)/lambdat1)*(mr_y1(NDtstar+NDtau)/lambdat2);
        ws(NDtstar) = 1/(lambdat1*lambdat2);
        %}
    end
        
    %S(NDtau+1) = sum(S_NDtau)/mean(ws); % TESTING
    S(NDtau+1)= sum(S_NDtau);
        
end

% For convenience only
NDtau = 0:maxtau;

end
