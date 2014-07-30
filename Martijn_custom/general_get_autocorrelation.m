function [S, Sraw, w, NDtau] = general_get_autocorrelation(y1,lambda)
% This function calculatates the autocorrelation for a vector with
% y-values. It can weigh certain parts of this vector more or less, and 
% will do so by weighing vector W. If w=ones(N_ttotal), this function will
% give the same output as xcov(y1,'coeff').
%
% Note that S(tau=0) = S(1) b/c matlab starts indexing at 1
%
% If w=0, then w will be set to ones(N)

% Some general stuff
% ===
% Total nr points
N=length(y1); 
% Importantly, use distances to mean ("mean removed y1").
mr_y1 = y1 - mean(y1); % 

% Set appropriate weights 
% ===
% Create set of empty weights, fill next
w=zeros(length(y1)+1); % +1 b/c indexing starts at 1 in matlab    
%(Note that sum(w(:,NDtau+1) should equal length(y1) for correct 
% weighing.)
for NDtau = 0:N
    for NDtstar = 1:(N-NDtau)
        % For no bias, everything weighs as one.
        if lambda==0
            w(NDtstar,NDtau+1)=1;
        else
            % Get weighing, Daan method
            if lambda(NDtstar) == 1
                w(NDtstar,NDtau+1) = 1/lambda(NDtstar+NDtau);
            else
                w(NDtstar,NDtau+1) = .75/lambda(NDtstar+NDtau);
            end
            % Weighing, Nature Kiviet:
            %{
            lambdat1 = lambda(NDtstar); % count at point 1
            lambdat2 = lambda(NDtstar+NDtau); % count at point 2
            w(NDtstar,NDtau+1)=1/(2*lambdat1)+1/(2*lambdat2); % weighing factor
            %}
        end
    end
end         

% Actual correlation
% ===
% Give me 
% S(tau) = 1/((t_total-tau)) sum_0^(t_total-tau) [ (y1(tstar)y2(tstar+tau)) dtstar ]
% R(tau) = S(tau)/S(0)
% In discrete form (tau -> N*Dtau = NDtau) etc 
S = []; 
for NDtau = 0:length(mr_y1)

    S_NDtau = [];
    %ws = []; % TEST
    
    % slide the window
    for NDtstar = 1:(length(mr_y1)-NDtau)
        
        % multiply two values
        S_NDtau(end+1) = mr_y1(NDtstar)*mr_y1(NDtstar+NDtau)*w(NDtstar,NDtau+1);
        
        % multiply two values TESTING
        %{
        if lambda==0, lambdat1=1; else lambdat1=lambda(NDtstar); end
        if lambda==0, lambdat2=1; else lambdat2=lambda(NDtstar+NDtau); end
        S_NDtau(end+1) = (mr_y1(NDtstar)/lambdat1)*(mr_y1(NDtstar+NDtau)/lambdat2);
        ws(NDtstar) = 1/(lambdat1*lambdat2);
        %}
    end
    
    % calculate average and record    
    S(NDtau+1) = sum(S_NDtau)/sum(w(:,NDtau+1)); % +1 b/c matlab starts indexing at 1
    %S(NDtau+1) = sum(S_NDtau)/mean(ws); % TESTING
    Sraw(NDtau+1)= sum(S_NDtau);
        
end

% For convenience only
NDtau = 0:length(mr_y1);

end
