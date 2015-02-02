function X = SSA_direct_function(s,prop,x0,T_Array)
%SSA_direct_function Calculates gillespie simulation.
%   SSA_direct_function gives the trajectory X of reactions 
%
% Written by Johannes Keegstra
% 2015/02

X=NaN(2,length(T_Array));

%run loop
t=0;
x=x0;
iii=1;
while t<(max(T_Array))
    
    w=prop(x);
    
    taumin=-log(rand(1,1))/sum(w); %the first reaction that occurs is a random variable.
    ra=rand(1,1)*sum(w);
    i=1;
    while sum(w(1:i))<ra
        i=i+1;
    end
    
    X(:,1)=x;
    
    %if no reaction occurs, proceed to next time step.
    if t>T_Array(iii);
    iii=iii+1;
    X(:,iii)=x;
    
    end
    
    
    
    
    
    t=t+taumin;
    if t>(max(T_Array));
          break
    else
           x=x+s(:,i); %update reactions.
    end
    
  
end



end

