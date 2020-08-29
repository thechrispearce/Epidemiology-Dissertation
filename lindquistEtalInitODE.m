function X0 = lindquistEtalInitODE(M,typeD,paramD,epsilon);

% initialize the vector X0 for use in Lindquist etal SIS-CM ED
% approximation

X0 = zeros((M+1)*(M+2),1);
%c = 0;

%pmf = @(x)poisspdf(x,lambda);
pmf = @(x)probsD(typeD,paramD,x,0);

%for s = 0:M
%    c = c + pmf(s);
%end

for s = 0:M
    %s
    %{
    X0(lindquistEtalStateTfm(M,1,s,0)) = epsilon*pmf(s)/c;
    X0(lindquistEtalStateTfm(M,0,s,0)) = (1-epsilon)*pmf(s)/c-epsilon*s*pmf(s)/c;
    if s>0
        X0(lindquistEtalStateTfm(M,0,s-1,1)) = epsilon*s*pmf(s)/c;
    end
    %}
    
    X0(lindquistEtalStateTfm(M,1,s,0)) = epsilon*pmf(s);
    X0(lindquistEtalStateTfm(M,0,s,0)) = (1-epsilon)*pmf(s)-epsilon*s*pmf(s);
    if s>0
        X0(lindquistEtalStateTfm(M,0,s-1,1)) = epsilon*s*pmf(s);
    end
    
    
    %X0
end