function [x_est] = ViterbiAlgorithm(N,s1,s2,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
% N : Number of Symbols
% s1 : Constant vector for x_n = 1
% s2 : Constant vector for x_n = -1
% r : Baseband equivalent signal
%
%------------------------------------------------------------
% OUTPUT :
% x_est : estimation of vector x
%
%------------------------------------------------------------
% Variable Naming Convention :
%
% 1 --> 3π/2 
% 2 --> π
% 3 --> π/2
% 4 --> 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializations
w12 = zeros(1,N); w32 = w12; w14 = w12; w34 = w12; w21 = w12; w41 = w12; w43 = w12; w1 = w12; w2 = w12; w3 = w12; w4 = w12; i1 = w12; i2 = w12; i3 = w12; i4 = w12; x_est = w12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward-Pass

i2(1) = -1;     % Set these purposely to -1 since on the very first step there is no transition from 0 -> 0 or from 0 -> π
i4(1) = -1;
w1(1) = real(r(:,1)'*s2);   % w2 from the trellis diagram
w3(1) = real(r(:,1)'*s1);   % w1 from the trellis diagram

for n=2:N
    if(mod(n,2)==0)
        w12(n) = real(r(:,n)'*s2*exp(1j*3*pi/2));
        w32(n) = real(r(:,n)'*s1*exp(1j*pi/2));
        w14(n) = real(r(:,n)'*s1*exp(1j*3*pi/2));
        w34(n) = real(r(:,n)'*s2*exp(1j*pi/2));
        
       x_est = zeros(1,N); t1 = w12(n)+w1(n-1);
        t2 = w32(n)+w3(n-1);
        [w2(n),i2(n)] = max([t1 0 t2 0]);
        
        t1 = w14(n)+w1(n-1);
        t2 = w34(n)+w3(n-1);
        [w4(n),i4(n)] = max([t1 0 t2 0]);
    else
        w21(n) = real(r(:,n)'*s1*exp(1j*pi));
        w41(n) = real(r(:,n)'*s2);
        w23(n) = real(r(:,n)'*s2*exp(1j*pi));
        w43(n) = real(r(:,n)'*s1);
        
        t1 = w21(n)+w2(n-1);
        t2 = w41(n)+w4(n-1);
        [w1(n),i1(n)] = max([0 t1 0 t2]);
        
        t1 = w23(n)+w2(n-1);
        t2 = w43(n)+w4(n-1);
        [w3(n),i3(n)] = max([0 t1 0 t2]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backward-Pass

% Keeping only one of the paths that gives the highest weight sum
if(mod(n,2)~=0)
    [~,path(N+1)] = max([w1(n) 0 w3(n) 0]);
else
    [~,path(N+1)] = max([0 w2(n) 0 w4(n)]);
end

path(1) = 4;
for n=N:-1:1
    if(mod(n,2)~=0 && n~=1)
        [~,i] = max([w1(n) 0 w3(n) 0]);
        in = [i1(n) 0 i3(n) 0];
        path(n) = in(i);
    elseif(mod(n,2)==0 && n~=1)
        [~,i] = max([0 w2(n) 0 w4(n)]);
        in = [0 i2(n) 0 i4(n)];
        path(n) = in(i);
    end
    
    if(path(n)-path(n+1)==-1)
        x_est(n) = -1;
    elseif(path(n)-path(n+1)==1)
        x_est(n) = 1;
    elseif(path(n)-path(n+1)==-3)
        x_est(n) = 1;
    elseif(path(n)-path(n+1)==3)
        x_est(n) = -1;
    end  
    
end

end