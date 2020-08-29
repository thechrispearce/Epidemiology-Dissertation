function X0 = PastorInitODE(M,typeD,paramD,epsilon)

% initialize the vector X0 for use in Pastoras
% approximation

X0 = zeros(M+1,1);
%c = 0;

%pmf = @(x)poisspdf(x,lambda);
pmf = @(x)probsD(typeD,paramD,x,0);

%for s = 0:M
%    c = c + pmf(s);
%end

for k = 0:M    
    X0(k+1) = epsilon;
end
    
    
    %X0
end