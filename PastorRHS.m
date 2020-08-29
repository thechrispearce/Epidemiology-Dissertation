function dy=PastorRHS(y,beta,gamma,M,typeD,paramD)

% calculate the right hand side of the Lindquist etal ODE model for the
% SIS-CM epidemic.

%first calculate the values of the big fraction
pmf = @(x)probsD(typeD,paramD,x,0);

meancomp=0;

for k=0:M
    meancomp=meancomp+k*pmf(k);
end

THETA=0;

for k=0:M
    THETA=THETA+(y(k+1)*k*pmf(k));
end

THETA=THETA/meancomp;

dy=zeros((M+1),1);

for k=0:M

    dy(k+1)=-gamma*y(k+1) + beta*k*(1-y(k+1))*THETA;
    
end